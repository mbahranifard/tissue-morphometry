close all hidden
clear all
clc

cd('C:\Users\mbahranifard3\Documents\Master_2019\Matlab\Stats Analysis\TM Cellularity\');
CellRaw = readtable ('ECMTurnOverUpdate.xlsx','Sheet',2, 'ReadVariableNames',false);
CellRaw = rows2vars(CellRaw);
CellRaw(:,1)=[];

%% temporary
celltemp  = CellRaw;
for i = 1:size(CellRaw,1)
    if ismember(CellRaw.Group(i),"1")
        celltemp{i,4}="Tg hAMSC Mid";
    elseif ismember(CellRaw.Group(i),"2")
        celltemp{i,4}="Tg Sham Mid";
    end
    if ismember(CellRaw.Annotator(i),"0")
        celltemp{i,5}="A";
    elseif ismember(CellRaw.Annotator(i),"1")
        celltemp{i,5}="B";
    end
end
writetable(celltemp,'C:\Users\mbahranifard3\Desktop\ECMturnover.csv');

%%

% CellRaw{:,6}=CellRaw.Var3./CellRaw.Var4;
% CellRaw.Properties.VariableNames={'Animal','Section','ECM','IW','Group','Ratio'};
CellRaw.Properties.VariableNames={'Animal','Section','Ratio','Group','Annotator'};
CellRaw.Section=transpose([1:size(CellRaw.Section,1)]);
G = groupsummary(CellRaw,{'Group','Animal'});
Gsub = groupsummary(G,{'Group'})
raw=CellRaw;


    figure('Units','pixels','WindowStyle','normal','Position',[50,50,700,700]);
    
%     raw = rmmissing(raw, 1);
    
    [secondarygrp,gpnum,IN2] = findgroups (raw.Group, raw.Animal);
%         splitknot =splitapply(@mean,raw.Dist,secondarygrp);
%         amatknot =splitapply(@std,raw.Dist,secondarygrp);
    [primarygrp,prgpnum] = findgroups (raw.Group);
    
    

        splitknot =splitapply(@mean,raw.Ratio,primarygrp);
        func = @(x) std(x)/sqrt(length(x));
        amatknot =splitapply(@std,raw.Ratio,primarygrp);
     grps=unique(secondarygrp);
    
     
%     for k = 1:size(grps)
%     trashvar = primarygrp;
%     trashvar (primarygrp == 1) =2;
%     trashvar (primarygrp == 2) =1;
%     primarygrp = trashvar;
    

%     mainx=dummyvar(categorical(prgpnum));

argvar = [2 1];

    bp=bar(argvar, splitknot,'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',2);
    
    colorvec = lines(numel(gpnum));
    hold all
    errorbar(argvar,splitknot,amatknot,'linestyle','none','color','k','linewidth',1.7)
%     end
%     secx=dummyvar(categorical(pri));
%     secx = secx *[1:11]';
    for i = 1:numel(gpnum)
    selection = find(secondarygrp == i);  
    if numel(unique(primarygrp(selection)))~=1
        sprintf("ridi")
    end
%     scatter(primarygrp(selection),raw{selection,4},70,'jitter','on','jitterAmount',0.25,'markerfacecolor',colorvec(i,:),'markeredgecolor','none')
    scatter(repelem(argvar(unique(primarygrp(selection))),numel(selection)),raw.Ratio(selection),120,'jitter','on','jitterAmount',0.35,'markerfacecolor',colorvec(i,:),'markeredgecolor','none','markerfacealpha',.75)

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
set(gca, 'xtick',[1:1:numel(argvar)],'xticklabel',{'Mid-term Sham' 'Mid-term hAMSC'}, 'XTickLabelRotation',45)
ylabel(['Basement Membrane Ratio'],'fontsize', 20)
ax.XAxis.FontSize =25;
ax.YAxis.FontSize =25;

%% fitLME
% %adding section ID
% uniqueElements = unique(CellRaw.Var5);
% numberingVector = zeros(size(CellRaw, 1), 1);
% for i = 1:numel(uniqueElements)
%     matchingRows = strcmp(CellRaw.Var5, uniqueElements{i});
%     numberingVector(matchingRows) = 1:sum(matchingRows);
% %     find(matchingRows)
% end
% CellRaw.Numbering = numberingVector;

%%
% clear all
% close all
% clc
% ds = readtable ('C:\Users\mbahranifard3\Desktop\Delivery method optimization\Delivery.csv');

% lme = fitlme(CellRaw,'ratio ~ group  + (1|group:eye)+ (group-1|group:eye)');
% lme = fitlme(CellRaw,'ratio ~ Var6  + (1|Var5)+ (Var6-1|Var5)');
CellRaw.Group = nominal(CellRaw.Group);
CellRaw.Animal = nominal(CellRaw.Animal);
CellRaw.Section = nominal(CellRaw.Section);
CellRaw.Annotator = nominal(CellRaw.Annotator);

% lme = fitlme(CellRaw,'ratio ~ Var6  + (Var6|Var5)+(Var6|Var5:Numbering)');
% formula='Ratio ~ Group  + (1|Animal)+(1|Animal:Section)';
formula='Ratio ~ Group  + (1|Animal)+(1|Animal:Section)+(1|Annotator)';
lme = fitlme(CellRaw,formula,'DummyVarCoding', 'reference');
%%test random effect
[~,~,STATS] = randomEffects(lme); % Compute the random-effects statistics (STATS)
STATS.Level = nominal(STATS.Level);
K = zeros(length(STATS),1);
K(STATS.Level == '10/9/2005') = 1;
pVal = coefTest(lme,[0 0 0 0 0 0 0 0 0],0,'REContrast',K')

% % %%%%
% lme = fitlme(CellRaw,'ratio ~ Var6  + (1|Var5)');
% lme = fitlme(CellRaw,'ratio ~ Var6');
% fitanov= anova(lme);
% comparisonResults = multcompare(fitanov);
num_fixed_effects = numel(unique(CellRaw.Group));
% Create an empty cell array to store the results of the t-tests and contrasts
all_results = zeros(num_fixed_effects);
[fixcoff,fixnames,stats] = fixedEffects(lme);
fixnames.Name(argvar)=fixnames.Name;
fixcoff(argvar)=fixcoff;
fixcoff(2)=0;
SElme = 1:10;
SElme(argvar) = stats.SE;
% Loop through each pair of fixed effects

%%
% very important about using the contrast vectors: in some cases the
% Y-intercept is already subtracted from the coefficient values. Make sure
% to add that back to coeffcieint before doing the following analysis
%%
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

