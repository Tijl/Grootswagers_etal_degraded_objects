%% prepare
subdir = '~/DATA/MEGBlurry/';
files = dir([subdir '*_200Hz.mat']);
nsubjects = 20;
for i=1:nsubjects
    subjects{i} = files(i).name(19:20);
end

% load all stimuli
imfiles = [];
for b={'Clear','Blurred'}
    for a={'Inanimate','Animate'}
        x=strsplit(ls(sprintf('~/PHD/Repository/Experiment1/StimuliMEG/%s/%s/*.png',b{1},a{1})));
        imfiles = [imfiles x(1:end-1)];
    end
end

load FQlevels

images = [];
for e=1:96 
    r = {};
    r.filepath = imfiles{e};
    r.image = imread(r.filepath);
    r.label = regexp(r.filepath,'/(\w+).png','tokens');
    r.label = r.label{1}{1};
    r.animate = ismember(e,[25:48 48+(25:48)]);
    r.blurred = e>48;
    r.number = e-r.blurred*48;
    r.exemplarnumber = e;
    r.orifilepath = regexp(r.filepath,'/(\w+/\w+.png)','tokens');
    r.orifilepath = sprintf('../Pictures/%s',r.orifilepath{1}{1});
    
    r.imf=fftshift(fft2(r.image));
    r.impf=abs(r.imf).^2; %power-frequency
    r.imsize=size(r.image);
    r.rotpf=rotavg(r.impf);
    
    [r.oriimage,~,r.orialpha] = imread(r.orifilepath);
    
    i=find(cellfun(@(x) strcmpi(x,sprintf('Pictures/%s.png',r.label)),values.imnames));
    r.oriorder = i;
    r.FQ = values.FQ(i);
    
    images = [images r];
end

%% FQ levels
ani=[images([images(:).blurred] & [images(:).animate]).FQ];
ina=[images([images(:).blurred] & ~[images(:).animate]).FQ];

fprintf('\nAnimate\t\tMean blur=%.2f (sd=%.2f)\n',mean(ani),std(ani))
fprintf('Inanimate\tMean blur=%.2f (sd=%.2f)\n',mean(ina),std(ina))


%% image stats
figure(1);clf;a=gca;hold on
hh=[];
for b=0:1
    for ani=1:-1:0
        sub=[images.animate]==ani & [images.blurred]==b;
        dat=[images(sub).rotpf]';
        h=plot(0:images(1).imsize/2,dat,'Color',a.ColorOrder(1+ani+2*b,:),'LineWidth',1)
        hh(end+1)=h(1);
    end
end
a.XScale='log';a.YScale='log';
legend(hh,{'animate','inanimate','animate degraded','inanimate degraded'})

%% LBA fit figure

lbadata_acc = readtable('figure_acc.csv');
lbadata_acc.y = 100-lbadata_acc.y;
lbadata_rt = readtable('figure_rt.csv');
lbadata_rt.y = lbadata_rt.y*1000;

f=figure(1);clf;f.Position=[0 100 900 300];f.PaperPositionMode='Auto';

leftm=.07;
rightm=.005;
width = (1-2*leftm-2*rightm)/2;
bottom = .1;
height = .84;

sd = lbadata_acc.y(1:4)-lbadata_acc.y(9:12);
x1 = [1 2];y1 = [3 4];y2 = [1 2];
shift = .0;

%accuracy animate
a=axes('Position',[leftm,bottom,width*.5,height]);a.FontSize=16;a.Box='on';hold on

h(1) = errorbar(x1-shift,lbadata_acc.y(y1),sd(y1),'.-','LineWidth',3,'MarkerSize',30,'Color',a.ColorOrder(1,:));
h(2) = plot(x1+shift,lbadata_acc.y(4+y1),'o:','LineWidth',3,'Color',a.ColorOrder(2,:));
ylabel('Accuracy (%)')
xlim([.5 2.5])
ylim([80 100])
a.XTick=1:2;a.YTick=[80:5:95];
a.XTickLabel={'Clear','Degraded'};a.LineWidth=1.5;
leg=legend([h(1) h(2)],{'Data','LBA'});leg.Orientation='horizontal';leg.Box='off';leg.Title.String='Animate';
leg.FontSize=16;leg.Position=[a.Position(1) bottom+.85*height .5*width bottom];

t=title('A');t.Position(1)=-.05;t.FontSize=35;t.HorizontalAlignment='left';t.VerticalAlignment='middle';

%accuracy inanimate
a=axes('Position',[leftm+.5*width,bottom,.5*width,height]);a.FontSize=16;a.Box='on';hold on
h(1) = errorbar(x1-shift,lbadata_acc.y(y2),sd(y2),'.-','LineWidth',3,'MarkerSize',30,'Color',a.ColorOrder(1,:));
h(2) = plot(x1+shift,lbadata_acc.y(4+y2),'o:','LineWidth',3,'Color',a.ColorOrder(2,:));
leg=legend([h(1) h(2)],{'Data','LBA'});leg.Orientation='horizontal';leg.Box='off';leg.Title.String='Inanimate';
leg.FontSize=16;leg.Position=[a.Position(1) bottom+.85*height .5*width bottom];
xlim([.5 2.5])
ylim([80 100])
a.XTick=1:2;a.YAxisLocation='Right';a.YTickLabel=[];a.LineWidth=1.5;
a.XTickLabel={'Clear','Degraded'};


%rt

sd = lbadata_rt.y(1:12)-lbadata_rt.y(25:36);
x1 = [1 1 1 2 2 2];y1 = 7:12;y2 = 1:6;

%rt animate
a=axes('Position',[2*leftm+rightm+width,bottom,width*.5,height]);a.FontSize=16;a.Box='on';hold on

for i = 1:3
    h(1) = errorbar(x1([i 3+i])-shift,lbadata_rt.y(y1([i 3+i])),sd(y1([i 3+i])),'.-','LineWidth',3,'MarkerSize',30,'Color',a.ColorOrder(1,:));
    h(2) = plot(x1([i 3+i])+shift,lbadata_rt.y(12+y1([i 3+i])),'o:','LineWidth',3,'Color',a.ColorOrder(2,:));
    
    %text
    xx = x1(3+i)+.1;
    yy = lbadata_rt.y(y1(3+i));
    tt=[10 50 90];
    t = text(xx,yy,sprintf('%i^{th}',tt(i)),'VerticalAlignment','Middle','FontSize',12);
    
end
leg=legend([h(1) h(2)],{'Data','LBA'});leg.Orientation='horizontal';leg.Box='off';leg.Title.String='Animate';
leg.FontSize=16;leg.Position=[a.Position(1) bottom+.85*height .5*width bottom];
ylabel('RT (ms)')
xlim([.5 2.5])
ylim([300 750])
a.XTick=1:2;a.YTick=300:100:700;
a.XTickLabel={'Clear','Degraded'};a.LineWidth=1.5;

t=title('B');t.Position(1)=-.1;t.FontSize=35;t.HorizontalAlignment='left';t.VerticalAlignment='middle';

%rt inanimate
a=axes('Position',[2*leftm+rightm+1.5*width,bottom,width*.5,height]);a.FontSize=16;a.Box='on';hold on

for i = 1:3
    h(1) = errorbar(x1([i 3+i])-shift,lbadata_rt.y(y2([i 3+i])),sd(y2([i 3+i])),'.-','LineWidth',3,'MarkerSize',30,'Color',a.ColorOrder(1,:));
    h(2) = plot(x1([i 3+i])+shift,lbadata_rt.y(12+y2([i 3+i])),'o:','LineWidth',3,'Color',a.ColorOrder(2,:));
    
    %text
    xx = x1(3+i)+.1;
    yy = lbadata_rt.y(y2(3+i));
    tt=[10 50 90];
    t = text(xx,yy,sprintf('%i^{th}',tt(i)),'VerticalAlignment','Middle','FontSize',12);
    
end
leg=legend([h(1) h(2)],{'Data','LBA'});leg.Orientation='horizontal';leg.Box='off';leg.Title.String='Inanimate';
leg.FontSize=16;leg.Position=[a.Position(1) bottom+.85*height .5*width bottom];
xlim([.5 2.5])
ylim([300 750])
a.XTick=1:2;a.YTick=300:100:700;
a.XTick=1:2;a.YAxisLocation='Right';a.YTickLabel=[];a.LineWidth=1.5;
a.XTickLabel={'Clear','Degraded'};

saveas(gcf,'Figures/lbafitstats','png')
saveas(gcf,'Figures/lbafitstats','fig')



%% load the LBA parameters and get distances
lbadata = readtable('parameters.csv');
conv_animate = [
    [21    8    9    2    6   14];
    [17   19   16    3   15   20];
    [10   24   23   12    5   13];
    [22    1   11   18    7    4];
];
conv_inanimate = [
    [    1   10   14   18   24   15];
    [   17   19    8   11    4   20];
    [    6    5   16    7    2   12];
    [   23   22   13    9   21    3];
];
RTresults = [];rng(1);
for s=1:nsubjects
    fprintf('Predicting RT for s %i/%i\n',s,nsubjects)

    [data,B] = loaddata(subjects{s},200);
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar+48*B.blurred,32);

    animatelabel = ismember(avlabels,[25:48 48+(25:48)]);
    blurredlabel = avlabels>48;

    res = timeseriesdecoding(avdata(:,:,:),animatelabel,'cvfolds',0,...
        'timevect',data.timevect,'verbose',2,'windowsize',5,'parallel',0);
    res.subject = subjects{s};
    res.exemplarlabel=1:96;
    res.animatelabel=animatelabel;
    res.blurredlabel=blurredlabel;

    % get RT
    res.rt = B.rt;
    res.animatetrials=B.animate;
    res.blurredtrials=B.blurred;
    res.zrt = zscore(B.rt);
    res.medianrt = zeros(96,1);
    res.medianzrt = zeros(96,1);
    for b=0:1
    for e=1:48
        res.medianrt(e+b*48) = median(res.rt(B.exemplar==e & B.blurred==b));
        res.medianzrt(e+b*48) = median(res.zrt(B.exemplar==e & B.blurred==b));
    end
    end
    
    % get acc
    res.acc = B.correct;
    res.accuracy = zeros(96,1);
    for b=0:1
    for e=1:48
        res.accuracy(e+b*48) = mean(res.acc(B.exemplar==e & B.blurred==b));
    end
    end
    
    % get lba drift
    clear_animate=[];clear_inanimate=[];blurry_animate=[];blurry_inanimate=[];
    sub = strcmp(lbadata.s,sprintf('subj%i',s)) & strcmp(lbadata.C,'TRUE');
    for x=1:4
        xx = sub & strcmp(lbadata.IW,sprintf('%i',x));
        for y=1:6
            yy = strcmp(lbadata.IB,sprintf('%i',y));
            clear_animate(conv_animate(x,y)) = lbadata.y(xx & strcmp(lbadata.D,'clear') & strcmp(lbadata.S,'animate') & yy);
            clear_inanimate(conv_inanimate(x,y)) = lbadata.y(xx & strcmp(lbadata.D,'clear') & strcmp(lbadata.S,'inanimate') & yy);
            blurry_animate(conv_animate(x,y)) = lbadata.y(xx & strcmp(lbadata.D,'blurry') & strcmp(lbadata.S,'animate') & yy);
            blurry_inanimate(conv_inanimate(x,y)) = lbadata.y(xx & strcmp(lbadata.D,'blurry') & strcmp(lbadata.S,'inanimate') & yy);
        end
    end
    res.lbarate = [clear_inanimate clear_animate blurry_animate blurry_inanimate]';
    RTresults = [RTresults res];
end
save RTresults.mat RTresults
    
%% add correlations
load RTresults.mat
RTresults1 = RTresults;
RTresults=[];
for s=1:nsubjects
    fprintf('Predicting RT for s %i/%i\n',s,nsubjects)
    [data,B] = loaddata(subjects{s},200);
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar+48*B.blurred,32);

    animatelabel = ismember(avlabels,[25:48 48+(25:48)]);
    blurredlabel = avlabels>48;
    res = RTresults1(s);
    %correlate
    ctype='Spearman';
    
    res.rtcall = corr(-res.medianzrt,res.targetposterior,'type',ctype);
    res.rtcani = corr(-res.medianzrt(animatelabel==1),res.targetposterior(animatelabel==1,:),'type',ctype);
    res.rtcina = corr(-res.medianzrt(animatelabel==0),res.targetposterior(animatelabel==0,:),'type',ctype);
    
    res.acccall = corr(res.accuracy,res.targetposterior,'type',ctype);
    res.acccani = corr(res.accuracy(animatelabel==1),res.targetposterior(animatelabel==1,:),'type',ctype);
    res.acccina = corr(res.accuracy(animatelabel==0),res.targetposterior(animatelabel==0,:),'type',ctype);
    
    res.lbacall = corr(res.lbarate,res.targetposterior,'type',ctype);
    res.lbacani = corr(res.lbarate(animatelabel==1),res.targetposterior(animatelabel==1,:),'type',ctype);
    res.lbacina = corr(res.lbarate(animatelabel==0),res.targetposterior(animatelabel==0,:),'type',ctype);
    
    res.clear_lbacall = corr(res.lbarate(blurredlabel==0),res.targetposterior(blurredlabel==0),'type',ctype);
    res.clear_lbacani = corr(res.lbarate(animatelabel==1 & blurredlabel==0),res.targetposterior(animatelabel==1 & blurredlabel==0,:),'type',ctype);
    res.clear_lbacina = corr(res.lbarate(animatelabel==0 & blurredlabel==0),res.targetposterior(animatelabel==0 & blurredlabel==0,:),'type',ctype);
    
    res.blurry_lbacall = corr(res.lbarate(blurredlabel==1),res.targetposterior(blurredlabel==1),'type',ctype);
    res.blurry_lbacani = corr(res.lbarate(animatelabel==1 & blurredlabel==1),res.targetposterior(animatelabel==1 & blurredlabel==1,:),'type',ctype);
    res.blurry_lbacina = corr(res.lbarate(animatelabel==0 & blurredlabel==1),res.targetposterior(animatelabel==0 & blurredlabel==1,:),'type',ctype);
   
    res.t=81;
    
    res.sumtargetposterior = sum(res.targetposterior(:,1:res.t),2);
    
    RTresults = [RTresults res];
end
save RTresults.mat RTresults

%% Summarize behaviour

f=figure(1);clf;a=gca;a.FontSize=16;f.PaperSize=[10 8];f.PaperPosition=[0 0 10 8];

D = {[RTresults.medianzrt],[RTresults.accuracy],[RTresults.lbarate]};
N = {'median zrt','accuracy','drift rate'};

for i=1:3
    subplot(3,1,i);a=gca;hold on;a.FontSize=16;
    d = D{i};

    dat=[reshape(d(RTresults(1).animatelabel==1 & RTresults(1).blurredlabel==0,:),[],1),...
        reshape(d(RTresults(1).animatelabel==0 & RTresults(1).blurredlabel==0,:),[],1),...
        reshape(d(RTresults(1).animatelabel==1 & RTresults(1).blurredlabel==1,:),[],1),...
        reshape(d(RTresults(1).animatelabel==0 & RTresults(1).blurredlabel==1,:),[],1),...
        ];
    
    [h,bc]=hist(dat,10);
    h = h/sum(h(:,1));m=mean(dat);
    h(:,3:4) = -h(:,3:4);
    p=plot(bc,h,'LineWidth',2);
    xlim(minmax(bc(:)'))
    plot(a.XLim,0*a.XLim,'k--')
    ylim(minmax(h(:)'))
    a.YTickLabel=abs(a.YTick);
    l=line([m;m],.25*[0 0 0 0;a.YLim(2) a.YLim(2) a.YLim(1) a.YLim(1)]);set(l,'LineWidth',4);
    if i==2
        legend(p,{'Animate clear','Inanimate clear','Animate blurry','Inanimate blurry'},'Location','EO')
    end
    xlabel(['exemplar ' N{i}])
    ylabel('proportion')
end
corr(squeeze(mean(cat(3,D{1:3}),2)))
saveas(gcf,'Figures/rtstats','png')
saveas(gcf,'Figures/rtstats','fig')

%% Summarize behaviour in line plot with anova

f=figure(1);clf;f.Position=[10 10 700 300];f.PaperPositionMode='Auto';

D = {1000*[RTresults.medianrt],100*[RTresults.accuracy],[RTresults.lbarate]};

N = {'Median RT (ms)','Mean accuracy (%)','Mean drift rate (A.U.)'};

for i=1:3
    leftm=.08;
    rightm=.01;
    width = (1-3*leftm-3*rightm)/3;
    bottom = .1;
    height = .85;
    switch(i)
        case 1;a=axes('Position',[leftm bottom width height]);
        case 2;a=axes('Position',[2*leftm+rightm+width bottom width height]);
        case 3;a=axes('Position',[3*leftm+2*rightm+2*width bottom width height]);
    end
    hold on;a.FontSize=16;
    d = D{i};

    dat=[mean(d(RTresults(1).animatelabel==1 & RTresults(1).blurredlabel==0,:));...
        mean(d(RTresults(1).animatelabel==0 & RTresults(1).blurredlabel==0,:));...
        mean(d(RTresults(1).animatelabel==1 & RTresults(1).blurredlabel==1,:));...
        mean(d(RTresults(1).animatelabel==0 & RTresults(1).blurredlabel==1,:));...
        ]';
    
    
    T = array2table(dat,'VariableNames',{'Y1','Y2','Y3','Y4'});
    within = table({'Clear';'Clear';'Degraded';'Degraded'},{'Animate';'Inanimate';'Animate';'Inanimate'},'VariableNames',{'degraded','animate'});
    rm = fitrm(T,'Y1-Y4~1','WithinDesign',within);
    ranova(rm,'WithinModel','degraded*animate')
    
    y = d(:);g={};
    g{1} = vertcat(RTresults.blurredlabel);
    g{2} = vertcat(RTresults.animatelabel);
    
    m=reshape(mean(dat),2,2);
    sd=reshape(standarderror(dat),2,2);
    
    e1=errorbar(.98*[1 2],m(:,1),sd(:,1),'linestyle','none','Marker','.','MarkerSize',30,'Color',a.ColorOrder(1,:),'LineWidth',3);
    e2=errorbar(1.02*[1 2],m(:,2),sd(:,2),standarderror(dat(:,[3 4])),'linestyle','none','Marker','.','MarkerSize',30,'Color',a.ColorOrder(3,:),'LineWidth',3);
    
    a.YLim=[.985 1.04].*[min(m(:)-sd(:)) max(m(:)+sd(:))];
    a.XLim=[.5 2.5];
    legend({'Clear','Degraded'},'location','NE','FontSize',15)
    a.XTick=[1 2];
    a.XTickLabel={'Animate','Inanimate'};a.XTickLabelRotation=0;
    ylabel(N{i})
    t=title(char(64+i));t.Units='Normalized';t.Position=[-0.25 .93];t.FontSize=30;

end

fprintf('MEDIAN RT: %.2f (sd: %.2f)\n',median(D{1}(:)),std(D{1}(:)))
fprintf('MEDIAN ACC: %.2f (sd: %.2f)\n',mean(D{2}(:)),std(D{2}(:)))
fprintf('MEDIAN LBA: %.2f (sd: %.2f)\n',mean(D{3}(:)),std(D{3}(:)))


saveas(gcf,'Figures/rtstatsbar','png')
saveas(gcf,'Figures/rtstatsbar','fig')


%% plot All

cc={'all','ani','ina'};close all
f=figure(1);clf;f.Position=[1 55 1000 600];f.PaperPositionMode='Auto';
for plotnr=1:3

    timevect=RTresults(1).timevect;
    eval(sprintf('A = vertcat(RTresults.lbac%s);',cc{plotnr}))
    eval(sprintf('B = vertcat(RTresults.rtc%s);',cc{plotnr}))
    eval(sprintf('C = vertcat(RTresults.accc%s);',cc{plotnr}))
    
    switch plotnr
        case 1
            subplot('Position',[.055 .07 .53 .88])
        case 2
            subplot('Position',[.645 .56 .34 .39])
        case 3
            subplot('Position',[.645 .07 .34 .39])
    end
    
    a=gca;a.FontSize=16;hold on;a.Position
    
    
    ha=shadedErrorBar(timevect,A,{@nanmean,@standarderror},{'color',a.ColorOrder(5,:),'linestyle','-','linewidth',3},1);
    hb=shadedErrorBar(timevect,B,{@nanmean,@standarderror},{'color',a.ColorOrder(2,:),'linestyle','-.','linewidth',3},1);
    hc=shadedErrorBar(timevect,C,{@nanmean,@standarderror},{'color',a.ColorOrder(4,:),'linestyle',':','linewidth',3},1);
    if plotnr==1
        pa=plotstarvect(timevect,A,0,-0.065,{'color',ha.mainLine.Color});
        pb=plotstarvect(timevect,B,0,-0.07,{'color',hb.mainLine.Color});
        pc=plotstarvect(timevect,C,0,-0.075,{'color',hc.mainLine.Color});
    else
        pa=plotstarvect(timevect,A,0,-0.065,{'color',ha.mainLine.Color});
        pb=plotstarvect(timevect,B,0,-0.075,{'color',hb.mainLine.Color});
        pc=plotstarvect(timevect,C,0,-0.085,{'color',hc.mainLine.Color});
    end
        
    
    plot(timevect,0*timevect,'k--')
        
    leg=legend([ha.mainLine,hb.mainLine,hc.mainLine],{'Drift rate','RT','Accuracy'},'Location','NE');
    if plotnr==1
        leg.FontSize=16;
    else
        leg.FontSize=16;
    end
        
    a.YColor='k';
    a.XLim=minmax(timevect);
    a.YLim=[-.1 .5];
    a.YTick=[-.1:.1:.45];a.Box='on';
    xlabel('time (ms)');
    ylabel('Spearman''s \rho');

    [m,ii]=max([mean(A);mean(B);mean(C)]');
    s=[standarderror(A);standarderror(B);standarderror(C)]';
    s=s(ii);

    pstim=patch([0 66 66 0],-.095+.005*[-1 -1 1 1],[.5 .5 .5],'EdgeColor','k');
    pn=patch([-90 60 60],[.1 .1 a.YLim(2)],'w','EdgeColor','w');
    b=[];colors = [ha.mainLine.Color;hb.mainLine.Color;hc.mainLine.Color];
    for i=1:length(m)
        p=patch(-100+40*i+18*[-1 1 1 -1],[.1 .1 m(i) m(i)],colors(i,:));
        ln=line(mean(p.XData([1:2;1:2]')),[m(i)-s(i),m(i)+s(i)]);ln.Color='k';ln.LineWidth=4;
        ln=line([p.XData(2),timevect(ii(i))],[m(i) m(i)]);ln.Color=p.FaceColor;ln.LineWidth=2;ln.LineStyle='--';
        b=[b p];
    end
    ln=line([-100 60 60],[.1 .1 a.YLim(2)]);ln.Color='k';ln.LineWidth=1;
   
    switch plotnr
        case 1;text(70,.5,{'','','','All exemplars'},'FontSize',16,'VerticalAlign','middle');
        case 2;text(70,.5,{'','','','Animate','exemplars'},'FontSize',16,'VerticalAlign','middle');
        case 3;text(70,.5,{'','','','Inanimate','exemplars'},'FontSize',16,'VerticalAlign','middle');
    end
    if plotnr==1
        title(char(64+plotnr),'units','data','Position',[-135 .5],'FontSize',35,'VerticalAlign','middle')
    else
        title(char(64+plotnr),'units','data','Position',[-150 .5],'FontSize',35,'VerticalAlign','middle')
    end
    if plotnr==3
        saveas(gcf,sprintf('Figures/predicting_behavior'),'png')
        saveas(gcf,sprintf('Figures/predicting_behavior'),'fig')
    end
end

%% correlate decoding with correlation
load decodingresults.mat

dectraj = mean(vertcat(decodingALL(:).balancedpcorr))';
[R,P]=corr(dectraj,[
mean(vertcat(RTresults.rtcall));
mean(vertcat(RTresults.rtcani));
mean(vertcat(RTresults.rtcina));
mean(vertcat(RTresults.lbacall));
mean(vertcat(RTresults.lbacani));
mean(vertcat(RTresults.lbacina));
mean(vertcat(RTresults.acccall));
mean(vertcat(RTresults.acccani));
mean(vertcat(RTresults.acccina));
]','type','Spearman')

R=reshape(R,3,3);
f=figure(1);clf;f.Position=[1,1,400,150];f.PaperPositionMode='auto';
imagesc(R',[0.15 .85]);
colormap parula
c=colorbar();c.FontSize=16;
c.Label.Rotation=0;c.Label.HorizontalAlignment='left';c.Label.VerticalAlignment='middle';
a=gca;a.FontSize=16;a.Position=[0.2    0.150    0.65    0.83];
a.XTick=[1 2 3];a.XTickLabel={'All','Animate','Inanimate'};
a.YTick=[1 2 3];a.YTickLabel={'Mean RT','Drift rate','Accuracy'};

P=reshape(P,3,3);
for p1=1:3
    for p2=1:3
        t=sprintf('%.2f',R(p1,p2));
        if 9*P(p1,p2)<0.001;t=[t '**'];elseif 9*P(p1,p2)<0.01;t=[t '*'];elseif 9*P(p1,p2)<0.05;t=[t '-'];end
        tx=text(p1-.3,p2,t);
        tx.FontSize=20;tx.HorizontalAlignment='left';tx.VerticalAlignment='middle';
        if R(p1,p2)<.3;tx.Color='w';end
    end
end
saveas(gcf,sprintf('Figures/correlating_correlations_with_decoding'),'png')
saveas(gcf,sprintf('Figures/correlating_correlations_with_decoding'),'fig')


%% correlate correlations responses
R=corr([
mean(vertcat(RTresults.rtcall));
mean(vertcat(RTresults.rtcani));
mean(vertcat(RTresults.rtcina));
mean(vertcat(RTresults.lbacall));
mean(vertcat(RTresults.lbacani));
mean(vertcat(RTresults.lbacina));
mean(vertcat(RTresults.acccall));
mean(vertcat(RTresults.acccani));
mean(vertcat(RTresults.acccina));
]','type','Spearman');
f=figure(1);clf;cla;f.PaperPosition=[0 0 10 10];
imagesc(flipud(R),[0.3 1]);axis square;a=gca;a.FontSize=16;colormap(a,colormap('bone'));a.TickLength=[0 0];
a.XTickLabel={'combined','animate','inanimate'};a.XAxis.FontSize=10;
a.YTickLabel=flipud(a.XTickLabel);a.YTickLabelRotation=90;a.YAxis.FontSize=10;
tx=text([2 5 8],[10 10 10],{'reaction time','drift rate','accuracy'},'FontSize',18,'HorizontalAlignment','Center','VerticalAlignment','Top');
ty=text([0 0 0],[2 5 8],fliplr({'reaction time','drift rate','accuracy'}),'Rotation',90,'FontSize',18,'HorizontalAlignment','Center','VerticalAlignment','Bottom');
lx=line(repmat([0.5 3.5 6.5 9.5],4,1),repmat([0;10],2,4));set(lx,'LineWidth',3,'Color','k')
ly=line(repmat([0;10],2,4),repmat([0.5 3.5 6.5 9.5],4,1));set(ly,'LineWidth',3,'Color','k')
ta=[];for i=1:9;for j=1:9;ta(i,j)=text(i,j,sprintf('%.2f',R(i,10-j)),'HorizontalAlignment','Center','Color',(R(i,10-j)<.5)*[1 1 1]);end;end;
saveas(gcf,'Figures/correlating_correlations','png')
saveas(gcf,'Figures/correlating_correlations','fig')

%% bar plot of important ones
f=figure(2);clf;cla;f.PaperPosition=[0 0 4 5];a=gca;hold on;a.Position=[.2 .1 .7 .8];

idx = [12;42;72]; %animate v inanimate
P=R(idx);


bfun = @(x) [
    corr(mean(vertcat(RTresults(x).lbacani))',mean(vertcat(RTresults(x).lbacina))','type','Spearman'),
    corr(mean(vertcat(RTresults(x).rtcani))',mean(vertcat(RTresults(x).rtcina))','type','Spearman'),
    corr(mean(vertcat(RTresults(x).acccani))',mean(vertcat(RTresults(x).acccina))','type','Spearman')];

P = bfun(1:20)';gcp();rng(1);
CI = bootci(10000,{bfun,1:20},'Options',struct('UseParallel',1));

colors=[5 4 2];
for i=1:3
    b=bar(i,P(i),.5);
    b.FaceColor=a.ColorOrder(colors(i),:);
end
e=errorbar(1:3,P,P-CI(1,:),CI(2,:)-P);e.LineStyle='none';e.LineWidth=2;e.Color='k';

a.FontSize=12;
a.XTick=1:3;a.XTickLabel={'Drift rate','RT','Accuracy'};
ylabel({'Correlation between animate & inanimate time series','(Spearman''s \rho)'})
ylim([-.6 1])
saveas(gcf,'Figures/correlating_correlations_bar','png')
saveas(gcf,'Figures/correlating_correlations_bar','fig')

%% scatter images ordered by distance

pnames = {'one minus median normalized rt','mean accuracy','mean drift rate'};
predictor=3;
showim=1;
timevect = RTresults(1).timevect;
close 1;f=figure(1);clf;f.Position=[0 0 1300 800];f.PaperPositionMode='auto';
subject = 1:20; 
for ani=1:-1:0
    if ani
        [~,t]=max(mean(vertcat(RTresults.lbacani)));t=63;
    else
        [~,t]=max(mean(vertcat(RTresults.lbacina)));t=63;
    end
    if ani
        a = axes('Position',[.04 .05 .45 .94])
    else
        a = axes('Position',[.54 .05 .45 .94])
    end
    a=gca;a.FontSize=16;hold on;


    xlabel(sprintf('distance to boundary (percentile) at %ims',timevect(t)))
    ylabel(sprintf('%s (percentile)',pnames{predictor}))
    xlim([-.1 1.05]);
    ylim([-.1 1.05]);
    animates = RTresults(1).animatelabel==ani;
    blurries = RTresults(1).blurredlabel(animates);
    switch predictor
        case 1;YY = -[RTresults.medianzrt];C=vertcat(RTresults.rtcall);
        case 2;YY = [RTresults.accuracy];C=vertcat(RTresults.acccall);
        case 3;YY = [RTresults.lbarate];C=vertcat(RTresults.lbacall);
    end
    
    Y=rankTransform(mean(YY(animates,subject),2),1);

    XX = cat(3,RTresults.targetposterior);
    X = rankTransform(mean(XX(animates,t,subject),3),1);


    h=line([X(blurries) X(~blurries)]',[Y(blurries) Y(~blurries)]');
    set(h,'LineWidth',1.5)
    c=sign(X(blurries)-X(~blurries))==sign(Y(blurries)-Y(~blurries));
    set(h(~c),'Color',a.ColorOrder(2,:))
    set(h(c),'Color',a.ColorOrder(5,:))

    ex = [1:24 49:72]+24*ani;
    w=.025;
    for i=1:length(ex)
        s=plot(X(i),Y(i),'o','MarkerFaceColor',[1 1 1],'MarkerSize',40,'LineWidth',3);
        s.MarkerEdgeColor = a.ColorOrder(2*images(ex(i)).blurred+1,:);

        image('CData',flip(images(ex(i)).oriimage,1),'XData',X(i)+[-w w],'YData',Y(i)+[-w w],'AlphaData',flip(images(ex(i)).orialpha,1));
    end

    %legend
    leg=legend( [plot(2,1,'o','Color',a.ColorOrder(1,:),'MarkerSize',10,'LineWidth',3),  ...
        plot(2,1,'o','Color',a.ColorOrder(3,:),'MarkerSize',10,'LineWidth',3), ...
        ],{'Clear','Degraded'} ,'Position',[a.Position(1)+.022 .875 .08 .05] );
    % add dists
    a1=axes('position',[a.Position(1) a.Position(2) a.Position(3) a.Position(4)/10]);hold on;
    xr = 0:.01:1;
    plot(xr,normpdf(xr,mean(X(~blurries)),std(X(~blurries))),'LineWidth',2,'Color',a.ColorOrder(1,:));
    plot(xr,normpdf(xr,mean(X(blurries)),std(X(blurries))),'LineWidth',2,'Color',a.ColorOrder(3,:));
    a1.Visible='off';
    a1.XLim=a.XLim;
    a2=axes('position',[a.Position(1) a.Position(2) a.Position(3)/10 a.Position(4)]);hold on;
    yr = 0:.01:1;
    plot(yr-.05,normpdf(yr,mean(Y(blurries)),std(Y(blurries))),'LineWidth',2,'Color',a.ColorOrder(1,:));
    plot(yr-.05,normpdf(yr,mean(Y(~blurries)),std(Y(~blurries))),'LineWidth',2,'Color',a.ColorOrder(3,:));
    a2.Visible='off';
    a2.View=[90 90];
    a2.XLim=a.XLim;

    switch ani
        case 1;title(a,'A       Animate','FontSize',25,'HorizontalAlignment','Left','Position',[-0.18,1,0])
        case 0;title(a,'B       Inanimate','FontSize',25,'HorizontalAlignment','Left','Position',[-0.18,1,0]);
    end
end
saveas(gcf,sprintf('Figures/Predicting %s images',pnames{predictor}),'png')
%saveas(gcf,sprintf('Figures/Predicting %s images',pnames{predictor}),'fig')

%% nice visualization of distance
t=63;
load RTresults

dat = cat(3,RTresults(:).targetposterior);
dat = squeeze(dat(:,t,:));
dat = mean(dat,2);
dat(ismember(1:96,[25:48 73:96]),:) = -dat(ismember(1:96,[25:48 73:96]),:);
rng(10000);
f=figure(1);clf;f.Position=[1 1 1000 1000];f.PaperPositionMode='auto';
a=gca;a.Position=[.05 .05 .94 .94];a.FontSize=20;a.Box='on';hold on
r=(repmat(1:(length(dat(:))/4),2,2)')*1.1;

xlim([-35 57])
ylim(xlim)
db=plot(.8*xlim,.8*xlim,'k-.','LineWidth',4);
%move some overlapping exemplars parallel to the boundary
fprintf('\n\n\n')
for i=1:48
    fprintf('%s\t%i\t%i\n',images(i).label,images(i).number,images(i).number+48)
end
%move clear up
for i=[34 40 20 24 10 10 9 9 32];r(i,1)=r(i,1)+1;end
%move clear down
for i=[27 27 47 21 7 2 5 6 1 17 17];r(i,1)=r(i,1)-1;end
%move blur up
for i=[22 22 23 23 8 8 5 7 1 43 15 13];r(i,2)=r(i,2)+1;end
%move blur down
for i=[28 46 46 20 24 26 30 33 32 11];r(i,2)=r(i,2)-1;end

%all ex
w=1.5; %width of images
im=[];s=[];
for b=0:1
    sub = 1:96>48 == b;
    d = mean(dat(sub,:),2);
    for i=1:48
        x=d(i)+r(i,1+b);
        y=-d(i)+r(i,1+b);
        s(b+1,i)=plot(x,y,'o','MarkerFaceColor',[1 1 1],'MarkerSize',40,'LineWidth',3);
        set(s(b+1,i),'MarkerEdgeColor',a.ColorOrder(2*b+1,:));
        im(b+1,i)=image('CData',flip(images(i).oriimage,1),'XData',x+[-w w],'YData',y+[-w w],'AlphaData',flip(images(i).orialpha,1));
        %debug: show numbers
        %text(x,y,sprintf('%i',i),'Rotation',45,'FontSize',20,'VerticalAlign','Middle','HorizontalAlign','Center')
    end
end

leg=legend( [plot(100,100,'o','Color',a.ColorOrder(1,:),'MarkerSize',10,'LineWidth',3),  ...
    plot(100,100,'o','Color',a.ColorOrder(3,:),'MarkerSize',10,'LineWidth',3), ...
    ],{'Clear objects','Degraded objects'} ,'Location','NE','FontSize',20);
tl = 28; toff = 0.5;
text(tl-toff,tl+toff,'animate','Rotation',45,'FontSize',20,'VerticalAlign','Bottom')
text(tl+toff,tl-toff,'inanimate','Rotation',45,'FontSize',20,'VerticalAlign','Top')
text(-20-toff,-20+toff,'decision boundary','Rotation',45,'FontSize',20,'VerticalAlign','Bottom')
a.XTick=[];a.YTick=[];
xlabel('Dimension 1')
ylabel('Dimension 2')
print(gcf,'-dpng','-r90','Figures/DistanceVisualization')
saveas(gcf,'Figures/DistanceVisualization','fig')

%% swarmplot-like visualization of distance, sorted by shift
t=63;
load RTresults

dat = cat(3,RTresults(:).targetposterior);
dat = squeeze(dat(:,t,:));
dat = mean(dat,2);
dat(ismember(1:96,[25:48 73:96]),:) = -dat(ismember(1:96,[25:48 73:96]),:);
rng(10000);
f=figure(1);clf;f.Position=[1 1 932 936];f.PaperPositionMode='auto';
a=gca;a.Units='pixels';a.Position=[30 35 900 900];a.FontSize=20;a.Box='on';hold on
a.LineWidth=2;
r=(repmat(1:(length(dat(:))/4),2,2)')*1.1;
%shift
sh = (mean(dat(1:48,:),2)-mean(dat(49:96,:),2)).*((1:48>24)'*2-1);
[~,idx]=sort(sh);
%all ex
w=0.02; %width of figures
im=[];s=[];

a.XLim=[0 1];
a.YLim=[0 1];

m1=-min(mean(dat(:,:),2))+2;
m2=max(m1+mean(dat(:,:),2))+2;

ln=line([m1/m2 m1/m2],a.YLim,'Color',a.ColorOrder(2,:),'LineWidth',2,'LineStyle','--');

for b=0:1
    sub = 1:96>48 == b;
    d = (m1+mean(dat(sub,:),2))./m2;
    for i=1:48
        x=d(idx(i));
        y=(1+2*i)/100;
        s(b+1,i)=plot(x,y,'o','MarkerFaceColor',[1 1 1],'MarkerSize',40,'LineWidth',3);
        set(s(b+1,i),'MarkerEdgeColor',a.ColorOrder(2*b+1,:));
        im(b+1,i)=image('CData',flip(images(idx(i)).oriimage,1),'XData',x+[-w w],'YData',y+[-w w],'AlphaData',flip(images(idx(i)).orialpha,1));
        %debug: show numbers
        %text(x,y,sprintf('%i',i),'Rotation',45,'FontSize',20,'VerticalAlign','Middle','HorizontalAlign','Center')
    end
end

leg=legend( [plot(-1,-1,'o','Color',a.ColorOrder(1,:),'MarkerSize',10,'LineWidth',3),  ...
    plot(-1,-1,'o','Color',a.ColorOrder(3,:),'MarkerSize',10,'LineWidth',3), ...
    ln,...
    ],{'Clear objects','Degraded objects','Decision boundary'},'Location','NW','FontSize',20);
a.XTick=[];a.YTick=[];
xlabel('Distance to classifier boundary')
ylabel('Exemplar')
print(gcf,'-dpng','-r300','Figures/DistanceVisualization2')
saveas(gcf,'Figures/DistanceVisualization2','fig')
