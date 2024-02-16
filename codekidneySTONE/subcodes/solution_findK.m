if algorithm == 1
  [Smax, Sid] = max(Sil);
  NCopt = NC(Sid)
  if Smax < 0.3
      R = length(NC);
      R = ceil(R/2):R;
      [Tmax, Q] = max(Silmin(R));
      Sid = R(Q);
      NCopt2 = NC(Sid);

  end
  disp([NC;Sil;Silmin]);
  if id == 14
    [TP,FP] = solution_positive(refseq_exon,refseq_intron,labelid,nsubset,Sid)
  end
 
elseif  algorithm == 2
  [Smax, Sid] = max(Sil);
  NCopt = NC(Sid);
  fprintf('\n## Clustering solution searched by Affinity Propagation:\n');
  fprintf('  Optimal number of clusters is %d, Silhouette = %g,\n',NCopt,Smax);
  if Smax < 0.3
      R = length(NC);
      R = ceil(R/2):R;
      [Tmax, Q] = max(Silmin(R));
      Sid = R(Q);
      NCopt2 = NC(Sid);
    
  end
 
  disp([NC;Sil;Silmin]);
  if id == 14
    [TP,FP] = solution_positive(refseq_exon,refseq_intron,labelid,nsubset,Sid);
   
  end
    
else
    NCs = unique(labels);
    NCopt = length(NCs);
    Sid = 1;
    [C, labels] = ind2cluster(labels);
    [NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCopt,...
      iend,simatrix,nrow,type,cut);

    if id == 14
      [TP,FP] = solution_positive(refseq_exon,refseq_intron,labels,labels(end),Sid);

    end
end