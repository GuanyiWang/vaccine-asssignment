%Computing the number of vaccinated units in each group
%author: Toru Kitagawa and Guanyi Wang

clear all

% DATA INPUT AND PREPARATION
load('../simulation_multi_fixden.mat')
load('../greedy_multi_fixden_para2_bgt1.mat')
filename_results = '../numbers15_greedy21.mat';
times = length(adj_1);
%initialize parameters
numa_21 = 0;
numb_21 = 0;
numa_22 = 0;
numb_22 = 0;
numa_23 = 0;
numb_23 = 0;
numa_31 = 0;
numb_31 = 0;
numa_32 = 0;
numb_32 = 0;
numa_33 = 0;
numb_33 = 0;
jj = 1;
aa = 1;
bb = 1;
cc = 1;
dd = 1;
ee = 1;

for ii = 1:times
    while jj <= length(V_result_21{ii})
    numa_21 = numa_21 + a_2(V_result_21{ii}(jj));
    numb_21 = numb_21 + b_2(V_result_21{ii}(jj));
    jj = jj+1;
    end
    
    while aa <= length(V_result_22{ii})
    numa_22 = numa_22 + a_2(V_result_22{ii}(aa));
    numb_22 = numb_22 + b_2(V_result_22{ii}(aa));
    aa = aa+1;
    end  
    
    while bb <= length(V_result_31{ii})
    numa_31 = numa_31 + a_3(V_result_31{ii}(bb));
    numb_31 = numb_31 + b_3(V_result_31{ii}(bb));
    bb = bb+1;
    end
    
    while cc <= length(V_result_32{ii})
    numa_32 = numa_32 + a_3(V_result_32{ii}(cc));
    numb_32 = numb_32 + b_3(V_result_32{ii}(cc));
    cc = cc+1;
    end
end


    while dd <= length(V_result_23)
    numa_23 = numa_23 + a_2(V_result_23(dd));
    numb_23 = numb_23 + b_2(V_result_23(dd));
    dd = dd+1;
    end
    
    while ee <= length(V_result_33)
    numa_33 = numa_33 + a_3(V_result_33(ee));
    numb_33 = numb_33 + b_3(V_result_33(ee));
    ee = ee+1;
    end
    
ratio_21 = numa_21/(numa_21+numb_21);
ratio_22 = numa_22/(numa_22+numb_22);
ratio_23 = numa_21/(numa_23+numb_23);
ratio_31 = numa_31/(numa_31+numb_31);
ratio_32 = numa_32/(numa_32+numb_32);
ratio_33 = numa_33/(numa_33+numb_33);

% Save coefficients and other variables for graphs
save(filename_results);
