

function [true_positive,false_positive] = solution_positive(refseq_exon,refseq_intron,labelid,classend,Sid)
ap_exon = (labelid(1:end-1,Sid)~=classend); 
true_positive=sum(ap_exon.*refseq_exon)/sum(refseq_exon)
false_positive=sum(ap_exon.*refseq_intron)/sum(refseq_intron)