
function [AGPSO1_cg_curve,PSO_cg_curve,IPSO_cg_curve]=Optcomparison(Max_iteration,SearchAgents_no)

Function_name='F8'; 
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Best_score1,Best_pos1,AGPSO1_cg_curve]= AGPSO1(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
[Best_score4,Best_pos4,PSO_cg_curve]   = PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
[Best_score5,Best_pos5,IPSO_cg_curve]= IPSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

display(['The best solution obtained by AGPSO1 is : ', num2str(Best_pos1)]);
display(['The best optimal value obtained by AGPSO1 is : ', num2str(Best_score1)]);

display(['The best solution obtained by PSO is : ', num2str(Best_pos4)]);
display(['The best optimal value obtained by SPSO is : ', num2str(Best_score4)]);

display(['The best solution obtained by IPSO is : ', num2str(Best_pos1)]);
display(['The best optimal value obtained by IPSO is : ', num2str(Best_score1)]);

        



