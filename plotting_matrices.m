


density_matrix_matlab = 'density_matrix_matlab.csv';

M = csvread(density_matrix_matlab)

figure 
surf(M);

shading interp
%FaceColor interp
ylabel('threshold', 'FontSize',12);

xlabel('internalisation rate','FontSize',12);

%yticks([0 10 20 30 40 50])
%yticklabels({'0.05','0.5','1','1.5','2','2.5'})


yticks([0 2 4 6 8])
yticklabels({'0.005','0.015','0.025','0.035','0.045'});


% for the fixed chemoattractant concentration
%xticks([0 10 20 30 40 50])
%xticklabels({'0.1','1','2','3','4','5'})

% for cell-induced model
%xticks([0 10 20 30 40 50])
%xticklabels({'10','100','200','300','400','500'})


xticks([0 2 4 6 8])
xticklabels({'30','90','150','210','270'});

title('Furthest distance travelled by cells')



set(gca,'FontSize',24)

  