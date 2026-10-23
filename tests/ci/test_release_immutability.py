"""Run real upload-step Bash against a local release-asset API substitute."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
FAKE_GH = r'''#!/usr/bin/env python3
import os
from pathlib import Path
import sys

args = sys.argv[1:]
if args[:3] != ["release", "upload", "v-test"]:
    raise SystemExit("unexpected release operation")
repo = args.index("--repo")
if args[repo + 1] != "example/rdistnn":
    raise SystemExit("unexpected repository")
assets = [Path(name) for name in args[3:repo]]
destination = Path(os.environ["LOCAL_RELEASE_ASSETS"])
if "--clobber" not in args and any((destination / asset.name).exists() for asset in assets):
    raise SystemExit("release asset already exists")
for asset in assets:
    (destination / asset.name).write_bytes(asset.read_bytes())
'''


@unittest.skipUnless(shutil.which("bash"), "bash is required")
class ReleaseImmutabilityTests(unittest.TestCase):
    def workflow_body(self, workflow, step_name):
        lines = (ROOT / ".github/workflows" / workflow).read_text().splitlines()
        step = lines.index("      - name: " + step_name)
        run = next(index for index in range(step + 1, len(lines))
                   if lines[index] == "        run: |")
        body = []
        for line in lines[run + 1:]:
            if line.strip() and not line.startswith("          "):
                break
            body.append(line[10:] if line.strip() else "")
        script = "\n".join(body).replace(
            "${{ github.event.release.tag_name }}", "v-test",
        ).replace("${{ github.repository }}", "example/rdistnn").replace(
            "${{ steps.makezip.outputs.zip_path }}", "_dist/rdistnn-source.zip",
        )
        self.assertNotIn("${{", script, "unresolved upload-step expression")
        return script

    def test_release_and_manual_assets_bypass_push_cooldown(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gh = root / "gh"
            gh.write_text('''#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys

with Path(os.environ["GH_CALLS"]).open("a") as calls:
    calls.write(" ".join(sys.argv[1:]) + "\\n")
route = sys.argv[2]
if route.endswith("/actions/runs/100"):
    result = {"workflow_id": 1, "head_branch": "master",
              "head_repository": {"full_name": "example/rdistnn"}}
elif "/actions/workflows/1/runs?" in route:
    result = {"workflow_runs": [{"id": 99, "event": os.environ["EVENT_NAME"],
              "head_branch": "master", "head_repository": {"full_name": "example/rdistnn"},
              "created_at": "2999-01-01T00:00:00Z"}]}
elif "/actions/runs/99/jobs?" in route:
    result = {"jobs": [{"name": "zip-and-upload", "started_at": "2999-01-01T00:00:00Z",
                        "status": "completed", "conclusion": "success"}]}
else:
    raise SystemExit("unexpected API call: " + route)
print(json.dumps(result))
''')
            gh.chmod(0o755)
            for event, reason, should_run in (
                ("release", "release_publication", "true"),
                ("workflow_dispatch", "manual_dispatch", "true"),
                ("push", "cooled_down", "false"),
            ):
                with self.subTest(event=event):
                    output, summary, calls = [root / (event + suffix)
                                              for suffix in (".output", ".summary", ".calls")]
                    env = os.environ.copy()
                    env.update(PATH=str(root) + os.pathsep + env.get("PATH", ""),
                               EVENT_NAME=event, WINDOW_HOURS="12", TARGET_JOB_NAME="zip-and-upload",
                               CURRENT_RUN_ID="100", REF_NAME="master", HEAD_REF="",
                               GITHUB_OUTPUT=str(output), GITHUB_STEP_SUMMARY=str(summary),
                               GH_CALLS=str(calls))
                    result = subprocess.run(
                        ["bash", "-c", self.workflow_body("package-snapshot.yml", "Enforce workflow cooldown")],
                        cwd=root, env=env, capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertIn("should_run=" + should_run + "\n", output.read_text())
                    self.assertIn("reason=" + reason + "\n", output.read_text())
                    self.assertEqual(calls.exists(), event == "push")

    def test_new_assets_upload_and_existing_assets_remain_unchanged(self):
        uploads = (
            ("plugin-build.yml", "Attach native plugins to release", "_release_assets",
             ("linux.plugin", "windows.plugin", "macos-arm.plugin", "macos-intel.plugin", "SHA256SUMS.txt")),
            ("package-snapshot.yml", "Attach validated source archive to release", "_dist",
             ("rdistnn-source.zip", "SHA256SUMS-source.txt")),
        )
        for workflow, step, source, filenames in uploads:
            for existing in (False, True):
                with self.subTest(workflow=workflow, existing=existing), tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    (root / "bin").mkdir()
                    gh = root / "bin/gh"
                    gh.write_text(FAKE_GH)
                    gh.chmod(0o755)
                    payload, release = root / source, root / "release"
                    payload.mkdir()
                    release.mkdir()
                    for name in filenames:
                        (payload / name).write_bytes(b"new payload")
                        if existing:
                            (release / name).write_bytes(b"immutable original")
                    env = os.environ.copy()
                    env.update(PATH=str(gh.parent) + os.pathsep + env.get("PATH", ""),
                               LOCAL_RELEASE_ASSETS=str(release))
                    result = subprocess.run(
                        ["bash", "-c", self.workflow_body(workflow, step)], cwd=root, env=env,
                        capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode == 0, not existing,
                                     result.stdout + result.stderr)
                    self.assertEqual({path.name for path in release.iterdir()}, set(filenames))
                    for path in release.iterdir():
                        self.assertEqual(path.read_bytes(),
                                         b"immutable original" if existing else b"new payload")


if __name__ == "__main__":
    unittest.main()
