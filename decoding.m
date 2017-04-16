%% prepare
subdir = '~/DATA/MEGBlurry/';
files = dir([subdir '*_200Hz.mat']);
nsubjects = 20;
for i=1:nsubjects
    subjects{i} = files(i).name(19:20);
end


%% all comparisons
decoding10foldC = [];
decoding10foldB = [];
TOTALTIME = [];
for s=1:nsubjects
    starttime = tic;
    fprintf('Decoding s %i/%i\n',s,nsubjects)

    [data,B] = loaddata(subjects{s},200);
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar+48*B.blurred,4);

    animatelabel = ismember(avlabels,[25:48 48+(25:48)]);
    blurredlabel = avlabels>48;
    avlabels(blurredlabel) = avlabels(blurredlabel)-48;
    
    %arguments for timeseriesdecoding
    parallel=1;
    windowsize=5;
    verbose=0;
    
    %clear
    subset = blurredlabel==0;
    res = timeseriesdecoding(avdata(subset,:,:),animatelabel(subset),...
        'timevect',data.timevect,'verbose',verbose,'windowsize',windowsize,'parallel',parallel);
    decoding10foldC = [decoding10foldC res];
    %blurry
    subset = blurredlabel==1;
    res = timeseriesdecoding(avdata(subset,:,:),animatelabel(subset),...
        'timevect',data.timevect,'verbose',verbose,'windowsize',windowsize,'parallel',parallel);
    decoding10foldB = [decoding10foldB res];
    
    % write out
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i - writing results\n',s,nsubjects);
    save decodingresults10fold.mat decoding10fold*
    
    TOTALTIME(s) = toc(starttime); %#ok<SAGROW>
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i ',s,nsubjects)
    fprintf('- TIME: %s ', datestr(TOTALTIME(s)*1/24/3600,'DD-HH:MM:SS'));
    fprintf('- TOTALTIME: %s ', datestr(sum(TOTALTIME)*1/24/3600,'DD-HH:MM:SS'));
    fprintf('- ETA: %s\n',datestr(mean(TOTALTIME(1:s))*(nsubjects-s)*1/24/3600,'DD-HH:MM:SS'))
end


%% all comparisons
decodingC2C = [];
decodingB2B = [];
decodingC2B = [];
decodingB2C = [];
decodingALL = [];
TOTALTIME = [];
for s=1:nsubjects
    starttime = tic;
    fprintf('Decoding s %i/%i\n',s,nsubjects)

    [data,B] = loaddata(subjects{s},200);
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar+48*B.blurred,4);

    animatelabel = ismember(avlabels,[25:48 48+(25:48)]);
    blurredlabel = avlabels>48;
    avlabels(blurredlabel) = avlabels(blurredlabel)-48;
    
    %arguments for timeseriesdecoding
    parallel=1;
    windowsize=5;
    verbose=0;
    
    %first decode all
    contrast='ALL';
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i - %s ',s,nsubjects,contrast);tic;
    res = timeseriesdecoding(avdata(:,:,:),animatelabel(:),'exemplarlabels',avlabels(:),...
        'timevect',data.timevect,'verbose',verbose,'windowsize',windowsize,'parallel',parallel);
    res.subject = subjects{s};
    res.type = contrast;
    decodingALL=[decodingALL res];
    fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
        
    %then do the others
    for train='CB'
        trainset = blurredlabel==strcmp(train,'B');
        for test='CB'
            contrast = [train '2' test];
            fprintf('%s  ',datestr(now))
            fprintf('subject %i/%i - %s ',s,nsubjects,contrast);tic;
            testset = blurredlabel==strcmp(test,'B');
            
            if strcmp(train,test) %normal decoding
                res = timeseriesdecoding(avdata(trainset,:,:),animatelabel(trainset),'exemplarlabels',avlabels(trainset),...
                    'timevect',data.timevect,'verbose',verbose,'windowsize',windowsize,'parallel',parallel);
            else %cross-training
                res = timeseriesdecoding(avdata(trainset,:,:),animatelabel(trainset),'exemplarlabels',avlabels(trainset),...
                    'testdata',avdata(testset,:,:),'testlabels',animatelabel(testset),'testexemplarlabels',avlabels(testset),...
                    'timevect',data.timevect,'verbose',verbose,'windowsize',windowsize,'parallel',parallel);
            end
            res.subject = subjects{s};
            res.type = contrast;
            res.trainset = trainset;
            res.testset = testset;
            %store
            eval(sprintf('decoding%s = [decoding%s res];',contrast,contrast))
            %print progress
            fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
        end
    end
    
    % write out
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i - writing results\n',s,nsubjects);
    save decodingresults.mat decoding*
    
    TOTALTIME(s) = toc(starttime); %#ok<SAGROW>
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i ',s,nsubjects)
    fprintf('- TIME: %s ', datestr(TOTALTIME(s)*1/24/3600,'DD-HH:MM:SS'));
    fprintf('- TOTALTIME: %s ', datestr(sum(TOTALTIME)*1/24/3600,'DD-HH:MM:SS'));
    fprintf('- ETA: %s\n',datestr(mean(TOTALTIME(1:s))*(nsubjects-s)*1/24/3600,'DD-HH:MM:SS'))
end

%%


%% plot the decoding results
load decodingresults.mat 
C2C = vertcat(decodingC2C.balancedpcorr);
B2B = vertcat(decodingB2B.balancedpcorr);
C2B = vertcat(decodingC2B.balancedpcorr);
B2C = vertcat(decodingB2C.balancedpcorr);
ALL = vertcat(decodingALL.balancedpcorr);
timevect = decodingC2C(1).timevect;

f=figure(1);clf;hold on;a=gca;a.FontSize=16;f.Position=[0 0 800 600];f.PaperPositionMode='Auto';
a.Position = [.06 .07 .91 .9];
ha=shadedErrorBar(timevect,C2C,{@mean,@standarderror},{'color',a.ColorOrder(1,:),'linestyle','-','linewidth',3},1);
hb=shadedErrorBar(timevect,B2B,{@mean,@standarderror},{'color',a.ColorOrder(3,:),'linestyle',':','linewidth',3},1);
hc=shadedErrorBar(timevect,ALL,{@mean,@standarderror},{'color',a.ColorOrder(2,:),'linestyle','-.','linewidth',3},1);
%hc=shadedErrorBar(timevect,C2B,{@mean,@standarderror},{'color',a.ColorOrder(3,:),'linewidth',2},1);
%hd=shadedErrorBar(timevect,B2C,{@mean,@standarderror},{'color',a.ColorOrder(4,:),'linewidth',2},1);
pa=plotstarvect(timevect,C2C,.5,0.45,{'color',ha.mainLine.Color});
pb=plotstarvect(timevect,B2B,.5,0.445,{'color',hb.mainLine.Color});
pc=plotstarvect(timevect,ALL,.5,0.44,{'color',hc.mainLine.Color});
%pc=plotstarvect(timevect,C2B,.5,0.44,{'color',hc.mainLine.Color});
%pd=plotstarvect(timevect,B2C,.5,0.435,{'color',hd.mainLine.Color});
pp=plotstarvect(timevect,C2C-B2B,0,0.435,{},'both');
plot(timevect,.5+0*timevect,'k--')
pstim=patch([0 66 66 0],.42+.01*[-1 -1 1 1],[.5 .5 .5],'EdgeColor','k');

leg=legend([ha.mainLine,hb.mainLine,hc.mainLine],...
    {'Clear','Degraded','Combined'},'Location','NW');
leg.FontSize=16;
xlabel('time (ms)')
ylabel('Decoding accuracy (% correct)')
ylim([.42 .75])
a.YTickLabel=100*a.YTick;
saveas(gcf,'Figures/decoding_animacy','png')
saveas(gcf,'Figures/decoding_animacy','fig')




%% plot the 10fold decoding results
load decodingresults.mat 
load decodingresults10fold.mat
C2C = vertcat(decodingC2C.balancedpcorr);
B2B = vertcat(decodingB2B.balancedpcorr);
C10 = vertcat(decoding10foldC.balancedpcorr);
B10 = vertcat(decoding10foldB.balancedpcorr);
timevect = decodingC2C(1).timevect;

f=figure(1);clf;hold on;a=gca;a.FontSize=16;f.PaperSize=[10 8];f.PaperPosition=[0 0 10 8];
% ha=shadedErrorBar(timevect,C2C,{@mean,@standarderror},{'color',a.ColorOrder(1,:),'linewidth',1,'linestyle','--'},1);
% hb=shadedErrorBar(timevect,B2B,{@mean,@standarderror},{'color',a.ColorOrder(2,:),'linewidth',1,'linestyle','--'},1);
hc=shadedErrorBar(timevect,C10,{@mean,@standarderror},{'color',a.ColorOrder(1,:),'linewidth',2},1);
hd=shadedErrorBar(timevect,B10,{@mean,@standarderror},{'color',a.ColorOrder(3,:),'linewidth',2},1);


%pa=plotstarvect(timevect,C2C,.5,0.45,{'color',ha.mainLine.Color});
%pb=plotstarvect(timevect,B2B,.5,0.445,{'color',hb.mainLine.Color});
pc=plotstarvect(timevect,C10,.5,0.44,{'color',hc.mainLine.Color});
pd=plotstarvect(timevect,B10,.5,0.435,{'color',hd.mainLine.Color});
pp=plotstarvect(timevect,C10-B10,0,0.43,{},'both');
plot(timevect,.5+0*timevect,'k--')

legend([hc.mainLine,hd.mainLine],...
    {'Clear','Degraded'},'Location','NW');
xlabel('time (ms)')
ylabel('Decoding accuracy (% correct)')
ylim([.4 .75])
a.YTickLabel=100*a.YTick;
saveas(gcf,'Figures/decoding_animacy_10fold','png')
saveas(gcf,'Figures/decoding_animacy_10fold','fig')

%% plot the difference
load decodingresults.mat 
load decodingresults10fold.mat
C2C = vertcat(decodingC2C.balancedpcorr);
B2B = vertcat(decodingB2B.balancedpcorr);
C10 = vertcat(decoding10foldC.balancedpcorr);
B10 = vertcat(decoding10foldB.balancedpcorr);
timevect = decodingC2C(1).timevect;

f=figure(1);clf;hold on;a=gca;a.FontSize=16;f.PaperSize=[10 8];f.PaperPosition=[0 0 10 8];
% ha=shadedErrorBar(timevect,C2C,{@mean,@standarderror},{'color',a.ColorOrder(1,:),'linewidth',2},1);
% hb=shadedErrorBar(timevect,B2B,{@mean,@standarderror},{'color',a.ColorOrder(2,:),'linewidth',2},1);
hc=shadedErrorBar(timevect,C10-C2C,{@mean,@standarderror},{'color',a.ColorOrder(1,:),'linewidth',2},1);
hd=shadedErrorBar(timevect,B10-B2B,{@mean,@standarderror},{'color',a.ColorOrder(2,:),'linewidth',2},1);
%pa=plotstarvect(timevect,C2C,.5,0.45,{'color',ha.mainLine.Color});
%pb=plotstarvect(timevect,B2B,.5,0.445,{'color',hb.mainLine.Color});
pc=plotstarvect2(timevect,C10,C2C,-0.014,{'color',hc.mainLine.Color},'both');
pd=plotstarvect2(timevect,B10,B2B,-0.016,{'color',hd.mainLine.Color},'both');
pe=plotstarvect2(timevect,C10-C2C,B10-B2B,-0.018,{'color','k'},'both');
plot(timevect,0+0*timevect,'k--')

legend([hc.mainLine,hd.mainLine],...
    {'Clear','Degraded'},'Location','NW');
xlabel('time (ms)')
ylabel('Difference between 10-fold and leave-one-exemplar-out')
ylim([-.02 .1])
saveas(gcf,'Figures/decoding_animacy_difference','png')
saveas(gcf,'Figures/decoding_animacy_difference','fig')





