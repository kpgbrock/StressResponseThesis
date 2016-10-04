load vhlUP;
orgseq = vhlUP.orgSeq(63:204);
[copt,V] = pholder(orgseq);
btorg = V*copt;
btorg = btorg/max(btorg);
aa = int2aa(1:20);
N = 1000;
results = zeros(N,1);
for i=1:N
    newseq = orgseq;
    newseq(randi(142,1)) = aa(randi(20,1));
    [copt,V] = pholder(newseq);
    btnew = V*copt;
    results(i) = sum(abs(btorg - btnew/max(btnew)));
end



    