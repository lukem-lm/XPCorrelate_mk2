%% PCDM2 - Oxygen edition
%Script with even more faffing about, largely following response to Angus'
%comments from the first dialogue of manuscripts (Mar2020)

%% Removing points with High O for visualisation

whereLowOX = datastack.EPMAO<3.5;
HlOx=H.*whereLowOX;
HlOx(HlOx==0)=NaN;
XPCcontourf(HlOx)

XPCcontourf(datastack.EPMAO)

nanmean(HlOx(:))

scatter(datastack.Phirefl(:)*180/pi,HlOx(:),'bx')
hold off 
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
%ylim([0 nanmean(H(:))+2*nanstd(H(:))])
xlim([0 90])
title('Declination angle against measured hardness with fit')
figname=['lowGBD Phi V H Figure FlippedFITCOS' ebsdname(1:(max(size(ebsdname)-4)))];
%% O vs H with colour based on Xcoord - original dataset

figure;
scatter3(datastack.EPMAO(:)*180/pi(),datastack.Hmatunsat(:),datastack.Phi(:));

temp=[datastack.EPMAO(:)*180/pi(),datastack.Hmatunsat(:),datastack.X(:)];
temp=sortrows(temp,3);
scatter(temp(:,1),temp(:,2),[],temp(:,3),'filled');
c=colorbar();
c.Label.String = 'Distance to GB /\mum';
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(H1(:))-5*nanstd(H1(:)) nanmean(H1(:))+5*nanstd(H1(:))])
%xlim([0 90])
title({'Declination angle against measured hardness',' with distance to grain boundary'})
figname=['Phi V H Figure_colours' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end

%% Normalising based on O

datastack.Hsmoothed=smoothdata(datastack.Hmat,1,'gaussian',4);
datastack.Hsmoothed=smoothdata(datastack.Hsmoothed,2,'gaussian',4);

%set up a new EPMA O signal which is smoothed, normalised, and reasonable.
xshift=mean(mean(datastack.EPMAOS(1:2,:)));

datastack.EPMAOSshift=datastack.EPMAOS-xshift;
datastack.EPMAOSshift(datastack.EPMAOSshift<0)=0;
datastack.EPMAOSshift(datastack.EPMAOSshift>20)=NaN;

%choose one model
Hpred= a(1)*((abs(datastack.EPMAOS-xshift)).^a(2));
%or overfit
Hpred=oxygenoverfit(a,datastack.EPMAOSshift);


XPCcontourf(datastack.BCebsd)

hold on
hist(datastack.EPMAOS(:)-xshift,1000)
XPCcontourf(datastack.EPMAOS)
min(datastack.EPMAOS(:))

HnormalisedO=datastack.Hsmoothed-Hpred;
XPCcontourf(HnormalisedO,'title',"Hardness, normalised by Oxygen",'cunits',"Hardness /GPa",'saveq',1)

meanx=mean(datastack.X,2);
meanHnormalisedO=nanmean(datastack.Hsmoothed,2);

scatter(datastack.X(:),datastack.Hsmoothed(:))
hold on
plot(meanx(:),meanHnormalisedO(:))

figure();
scatter(datastack.Phi(:)*180/pi,HnormalisedO(:))
xlim([0 90])
xlabel('Declination Angle /^{o}')
ylabel('Normalised Hardness /GPa')
%ylim([nanmean(H1(:))-5*nanstd(H1(:)) nanmean(H1(:))+5*nanstd(H1(:))])
title({'Declination angle against normalised (by O) measured hardness'})
figname=['Phi V H norm figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end
%% Pulling out invidual grains of interest

%first lets get rid of bakelite
setofinterest = datastack.H < 3.5 & datastack.Phi < 20*pi()/180;
setofinterest = datastack.H < 3.2 & datastack.Phi < 30*pi()/180 & datastack.Phi > 20*pi()/180;

setofinterest =  datastack.Phi2 < 8;
setofinterest =  datastack.GBDseccorr > 7;


Hofinterest=datastack.H.*setofinterest;
Hofinterest(Hofinterest==0)=NaN;

figure;
hplot=contourf(X,Y,Hofinterest,45,'LineColor','None');
title('Set of interest Hardness Nanoindentation Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['Hofinterest Figure paper' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')

%% Normalising H and plotting w.r.t GBD
Hpred=a(1)*cos(2*datastack.Phireflseccorr)+a(2);

datastack.Hsmoothed=smoothdata(datastack.H,1,'gaussian',2);
datastack.Hsmoothed=smoothdata(datastack.Hsmoothed,2,'gaussian',2);


Hnormalisedphi=datastack.Hsmoothed-Hpred;
XPCcontourf(Hpred,'title',"Predicted Hardness - smoothed",'cunits',"Hardness /GPa",'saveq',1)

XPCcontourf(Hnormalisedphi,'title',"Normalised Hardness",'cunits',"Hardness /GPa",'saveq',1)


%{
gmodsqrbsd=sqrt(gmo(:,:,1).^2+gmo(:,:,2).^2+gmo(:,:,3).^2)./GBD;
setofinterest =  datastack.GBD < 8;
figure;
XPCcontourf(gmodsqrbsd2)
gmodsqrbsd2=gmodsqrbsd.*setofinterest;
gmodsqrbsd2(gmodsqrbsd2==0)=NaN;
figure;
scatter(gmodsqrbsd2(:),Hnormalisedphi(:),[],datastack.GBD(:));
xlabel('Distance to GB /\mum')
ylabel('Normalised Hardness /GPa')
%}

figure;
scatter(datastack.GBDseccorr(:),Hnormalisedphi(:));
xlabel('Distance to GB /\mum')
ylabel('Normalised Hardness /GPa')
%ylim([nanmean(H1(:))-5*nanstd(H1(:)) nanmean(H1(:))+5*nanstd(H1(:))])
%xlim([0 90])
ylim([-2 2])
title({'GB distance against normalised measured hardness'})
figname=['GBD V H norm figure blankseccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end


%Binning the GBD and averaging. 

%make an average for all 10deg bins:
hgbdbinno=50;
gbdave=zeros(hgbdbinno,1);
Hgbdave=zeros(hgbdbinno,1);
stdevHgbdave=zeros(hgbdbinno,1);
countnum=zeros(hgbdbinno,1);
for i = 1:hgbdbinno
    wherei = datastack.GBDseccorr >(i-1)*max(datastack.GBDseccorr(:))/hgbdbinno & datastack.GBDseccorr <i*max(datastack.GBDseccorr(:))/hgbdbinno;
    gbdave(i)=(i-0.5)*max(datastack.GBDseccorr(:))/hgbdbinno;
    Hgbdave(i) = nanmean(Hnormalisedphi(wherei));
    stdevHgbdave(i) =nanstd(Hnormalisedphi(wherei));
    countnum(i)=nnz(wherei);
end
scatter(gbdave(:),Hgbdave(:))
hold on
errorbar(gbdave(:),Hgbdave(:),stdevHgbdave(:),'LineStyle','none');
xlabel('Distance to GB /\mum')
ylabel('Hardness /GPa')
%text(gbdave(:)+3,Hgbdave(:)+0.05,num2cell(countnum))
title('Hardness near GB, Averaged')
figname=['HnearGB Binned seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
hold off


%%
% Another way to visualise this: 
imgcomp=imshowpair(rot90(datastack.Phipred),rot90(datastack.Phi))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
figname=['Phipred vs Phi' ebsdname(1:(max(size(ebsdname)-4))) '.png'];
imwrite(imgcomp.CData,fullfile(resultsdir, figname))

imgcomp=imshowpair(rot90(datastack.Phipred),rot90(datastack.Phi),'diff')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
figname=['Phipred vs Phi2' ebsdname(1:(max(size(ebsdname)-4))) '.png'];
imwrite(imgcomp.CData,fullfile(resultsdir, figname))

%Now doing this with the second correction:
imgcomp=imshowpair(rot90(datastack.Phipredseccorr),rot90(datastack.Phiseccorr),'diff')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
figname=['Phipred vs Phi2seccorr' ebsdname(1:(max(size(ebsdname)-4))) '.png'];
imwrite(imgcomp.CData,fullfile(resultsdir, figname))


%% Correction trials - 
% This code is just to carry out a bunch of tests with 4, 5, 6 control
% points, and see how many points might have different assignment of
% grainId. 
gIDtrials.t1_4points=datastack.gID;
gIDtrials4p.rsquared;
gIDtrials5p.rsquared;
gIDtrials6p.rsquared=R2;
gIDtrialsall.gIDtrials4p=gIDtrials4p;
gIDtrialsall.gIDtrials5p=gIDtrials5p;
gIDtrialsall.gIDtrials6p=gIDtrials6p;

fieldloop=[{'t1'},{'t2'},{'t3'},{'t4'}];
for i=2:size(fieldloop,2)
    [datastack]=f_fixEBSDdistortion(ebsd,X,Y,fullres,microscope,primphase, EBSDREFq,resultsdir);
    datastack.Phi(datastack.Phi>(pi/2))=pi-datastack.Phi(datastack.Phi>(pi/2));
    f_writeEBSDdata(fullfile(resultsdir, [ebsdname(1:(max(size(ebsdname)-4))) 'corrected' currdate]),datastack)
    datastack=f_dist2grainb(resultsdir, ebsdname,datastack,currdate,ebsd);
    gIDtrials6p.(fieldloop{i})=datastack.gID;
end



fields = fieldnames(gIDtrials6p);
for k=1:numel(fields)
    gIDoftrials6(:,k)=gIDtrials6p.(fields{k})(:);
end
gIDtrialling(:,1)=gID7apr(:);
gIDtrialling(:,2)=datastack.gID(:);
gIDtrialling(:,3)=datastack.gIDseccorr(:);


% NOTES:https://uk.mathworks.com/matlabcentral/answers/406619-3d-coordinates-line-of-fit
X_ave=mean(gIDtrialling,1);            % mean; line of best fit will pass through this point  
dX=bsxfun(@minus,gIDtrialling,X_ave);  % residuals
C=(dX'*dX)/(size(gIDtrialling,1)-1);           % variance-covariance matrix of X
[R,D]=svd(C,0);             % singular value decomposition of C; C=R*D*R'
D=diag(D);
R2=D(1)/sum(D);
disp(R2)
plot3(gIDtrialling(:,1),gIDtrialling(:,2),gIDtrialling(:,3),'.k','MarkerSize',13)           % simulated noisy data
xlabel('gID from 6 points')
ylabel('gID from 9 points')
zlabel('gID from 9 points, after second correction')
% IF YOU HAVE 3 VARIABLES, YOU CAN DO THE FOLLOWING:
%{
% Visualize X and line of best fit
% -------------------------------------------------------------------------
% End-points of a best-fit line (segment); used for visualization only 
x=dX*R(:,1);    % project residuals on R(:,1) 
x_min=min(x);
x_max=max(x);
dx=x_max-x_min;
Xa=(x_min-0.05*dx)*R(:,1)' + X_ave;
Xb=(x_max+0.05*dx)*R(:,1)' + X_ave;
X_end=[Xa;Xb];
figure('color','w')
axis equal 
hold on
plot3(X_end(:,1),X_end(:,2),X_end(:,3),'-r','LineWidth',3) % best fit line 
set(get(gca,'Title'),'String',sprintf('R^2 = %.3f',R2),'FontSize',25,'FontWeight','normal')
view([20 20])
drawnow
%}


%%
figure;
plot(ebsdCorrectedLGG('Titanium'),ebsdCorrectedLGG('Titanium').orientations)
figure();

XPCcontourf(log(abs(Hnormalisedphi./datastack.GBDseccorr)))
XPCcontourf(mPrimeave./datastack.GBDseccorr)
XPCcontourf(datastack.GBD./datastack.GBDseccorr)

XPCcontourf(datastack.GBD)
XPCcontourf(datastack.GBDseccorr)

wherelowGBD=datastack.GBDseccorr<5;
H1=Hnormalisedphi.*wherelowGBD;
XPCcontourf(H1)


scatter(mPrimeave(:)./datastack.GBDseccorr(:),log(abs(H1(:)./datastack.GBDseccorr(:))))