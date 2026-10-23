version 17
// Independent finite-sample references: enumerate s-subsets and assign each
// selected tie its share among the subset's exact-distance boundary members.
// Production ranks, weight tables, block and deletion helpers are not used.
foreach fn in subset_weights subset_deletions distance_table distance_deletions {
    capture mata: mata drop rdistnn_test_`fn'()
}
mata:
real matrix rdistnn_test_subset_weights(real matrix DX,real scalar k,real scalar s)
{
    real scalar n,m,mask,j,cut,below,tied,total
    real colvector selected,dist,inside,equal
    real matrix W
    n=rows(DX); m=cols(DX); W=J(n,m,0); total=0
    for(mask=0;mask<2^n;mask++) {
        selected=mod(floor(mask:/(2:^(0::(n-1)))),2)
        if(sum(selected)!=s) continue
        total++
        for(j=1;j<=m;j++) {
            dist=sort(select(DX[,j],selected),1); cut=dist[k]
            inside=selected:&(DX[,j]:<cut); equal=selected:&(DX[,j]:==cut)
            below=sum(inside); tied=sum(equal)
            W[,j]=W[,j]+inside:/k+equal:*((k-below)/(tied*k))
        }
    }
    return(W:/total)
}
real matrix rdistnn_test_subset_deletions(real matrix DX,real matrix S,real scalar k,real scalar s)
{
    real scalar i,j,n,q
    real colvector keep
    real matrix W,L
    n=rows(DX); q=cols(S); L=J(n,cols(DX)*q,.)
    for(i=1;i<=n;i++) {
        keep=selectindex((1::n):!=i)
        W=rdistnn_test_subset_weights(DX[keep,],k,s)
        for(j=1;j<=cols(DX);j++) L[i,((j-1)*q+1)..(j*q)]=(S[keep,]'*W[,j])'
    }
    return(L)
}
real matrix rdistnn_test_distance_table(real matrix DX,real matrix R,real colvector table)
{
    real scalar i,j,n
    real colvector equal
    real matrix W
    n=rows(DX); W=J(n,cols(DX),.)
    for(j=1;j<=cols(DX);j++) for(i=1;i<=n;i++) {
        equal=(DX[,j]:==DX[i,j])
        W[i,j]=quadsum(table[select(R[,j],equal)])/sum(equal)
    }
    return(W)
}
real matrix rdistnn_test_distance_deletions(real matrix DX,real matrix R,real matrix S,real colvector table)
{
    real scalar i,j,n,q
    real colvector keep
    real matrix RR,W,L
    n=rows(DX); q=cols(S); L=J(n,cols(DX)*q,.)
    for(i=1;i<=n;i++) {
        keep=selectindex((1::n):!=i)
        RR=R[keep,]:-(R[keep,]:>J(n-1,1,1)*R[i,])
        W=rdistnn_test_distance_table(DX[keep,],RR,table)
        for(j=1;j<=cols(DX);j++) L[i,((j-1)*q+1)..(j*q)]=(S[keep,]'*W[,j])'
    }
    return(L)
}
end
