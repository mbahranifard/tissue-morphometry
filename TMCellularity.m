close all hidden
clear all
clc

cd('C:\Users\mbahranifard3\Documents\Master_2019\Matlab\Stats Analysis\TM Cellularity\');
CellRaw = readtable ('TMCellularity.xlsx','Sheet',2,'Format','auto');
G = groupsummary(CellRaw,{'Var6','Var5'});
Gsub = groupsummary(G,{'Var6'})
raw=CellRaw;


    figure('Units','pixels','WindowStyle','normal','Position',[50,50,700,700]);
    
    raw = rmmissing(raw, 1);
    
    [secondarygrp,gpnum,IN2] = findgroups (raw.Var5, raw.Var6);
%         splitknot =splitapply(@mean,raw.Dist,secondarygrp);
%         amatknot =splitapply(@std,raw.Dist,secondarygrp);
    [primarygrp,prgpnum] = findgroups (raw.Var6);
    
    

        splitknot =splitapply(@mean,raw(:,4),primarygrp);
        func = @(x) std(x)/sqrt(length(x));
        amatknot =splitapply(@std,raw(:,4),primarygrp);
     grps=unique(secondarygrp);
    
     
%     for k = 1:size(grps)
%     trashvar = primarygrp;
%     trashvar (primarygrp == 1) =2;
%     trashvar (primarygrp == 2) =1;
%     primarygrp = trashvar;
    

%     mainx=dummyvar(categorical(prgpnum));

argvar = [2 10 6 4 9 5 3 8 7 1];

    bp=bar(argvar, splitknot,'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',2);
    splitknot
    colorvec = lines(numel(gpnum));
    hold all
    errorbar(argvar,splitknot,amatknot,'linestyle','none','color','k','linewidth',1.2)
%     end
%     secx=dummyvar(categorical(pri));
%     secx = secx *[1:11]';
    for i = 1:numel(gpnum)
    selection = find(secondarygrp == i);  
    if numel(unique(primarygrp(selection)))~=1
        sprintf("ridi")
    end
%     scatter(primarygrp(selection),raw{selection,4},70,'jitter','on','jitterAmount',0.25,'markerfacecolor',colorvec(i,:),'markeredgecolor','none')
    scatter(repelem(argvar(unique(primarygrp(selection))),numel(selection)),raw{selection,4},70,'jitter','on','jitterAmount',0.25,'markerfacecolor',colorvec(i,:),'markeredgecolor','none')

    end  
    %% reporting stast
            oavg(argvar)=splitknot;
            ostd(argvar)=amatknot;
            numdif(argvar) = Gsub.GroupCount;
            CIcell = ostd./sqrt(numdif).*tinv(0.975,numdif-1); 
            Ystat = transpose([oavg;oavg-CIcell;oavg+CIcell])

    %%
    
ylim ([0 1])
ax = gca;
ax.FontSize = 18;
set(gca, 'xtick',[1:1:numel(argvar)],'xticklabel',{'WT' 'Het' 'ShortTg' 'ShortSC' 'MidTg' 'MidSC','LongTg','LongSc','ShortiPSC','MidiPSC'}, 'XTickLabelRotation',45)
ylabel(['Normalized TM Cell Count (cell/\mum)'],'fontsize', 20)
ax.XAxis.FontSize =25;
ax.YAxis.FontSize =25;

%% fitLME
%adding section ID
uniqueElements = unique(CellRaw.Var5);
numberingVector = zeros(size(CellRaw, 1), 1);
for i = 1:numel(uniqueElements)
    matchingRows = strcmp(CellRaw.Var5, uniqueElements{i});
    numberingVector(matchingRows) = 1:sum(matchingRows);
%     find(matchingRows)
end
CellRaw.Numbering = numberingVector;

%%
% clear all
% close all
% clc
% ds = readtable ('C:\Users\mbahranifard3\Desktop\Delivery method optimization\Delivery.csv');

% lme = fitlme(CellRaw,'ratio ~ group  + (1|group:eye)+ (group-1|group:eye)');
% lme = fitlme(CellRaw,'ratio ~ Var6  + (1|Var5)+ (Var6-1|Var5)');
CellRaw.Var6 = nominal(CellRaw.Var6);
CellRaw.Numbering = nominal(CellRaw.Numbering);
CellRaw.Var5 = nominal(CellRaw.Var5);

% lme = fitlme(CellRaw,'ratio ~ Var6  + (Var6|Var5)+(Var6|Var5:Numbering)');
formula='ratio ~ Var6  + (1|Var5)+(1|Var5:Numbering)';
lme = fitlme(CellRaw,formula,'DummyVarCoding', 'reference');
% lme = fitlme(CellRaw,'ratio ~ Var6  + (1|Var5)');
% lme = fitlme(CellRaw,'ratio ~ Var6');
% fitanov= anova(lme);
% comparisonResults = multcompare(fitanov);
num_fixed_effects = numel(unique(CellRaw.Var6));
% Create an empty cell array to store the results of the t-tests and contrasts
all_results = zeros(num_fixed_effects);
[fixcoff,fixnames,stats] = fixedEffects(lme);
fixnames.Name(argvar)=fixnames.Name;
fixcoff(argvar)=fixcoff;
fixcoff(2)=0;
SElme = 1:10;
SElme(argvar) = stats.SE;
% Loop through each pair of fixed effects
for i = 1:num_fixed_effects
    for j = i+1:num_fixed_effects
        
        % Define the contrast vector for this pair of fixed effects
        contrast_vector = zeros(1, num_fixed_effects);
        contrast_vector(i) = 1;  % Set weight 1 for the first fixed effect
        contrast_vector(j) = -1; % Set weight -1 for the second fixed effect
        % Create a contrast matrix
        contrast_matrix = contrast_vector;
        SEcomp = sqrt(SElme(i)^2+SElme(j)^2);
        t_stat = (contrast_matrix * fixcoff)/SEcomp;
        p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));
        % Store the results in the cell array
        all_results(i, j) = p_value;
    end
end
lmecomptbl = table(fixnames.Name,all_results)


F = fitted(lme);
R = residuals(lme);

%plot fitted line
[ypred,yCI,~] = predict(lme,CellRaw);

figure()
h1 = line((CellRaw.Var6),ypred);
hold on;
h2 = plot(CellRaw.Var6,yCI,'g-.');
gscatter(CellRaw.Var6,CellRaw.ratio,CellRaw.Var5)
hold on
gscatter(CellRaw.Var6,ypred,CellRaw.Var5,[],'o+x')
% plot(F,R,'bx')
% xlabel('Fitted Values')
% ylabel('Residuals')
% Apply Bonferroni correction manually

