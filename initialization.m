%_________________________________________________________________________%
%  Fungal Growth Optimizer source code demo 1.0               %
%  Developed in MATLAB R2019A                                      %
%                                                                         %
%                                                                         %
%   Main paper: Fungal Growth Optimizer: A Novel Nature-inspired Metaheuristic Algorithm for Stochastic Optimization                       %
%                  %
%                                                                         %
%_________________________________________________________________________%

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= length(ub); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
     Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;      
    end
end