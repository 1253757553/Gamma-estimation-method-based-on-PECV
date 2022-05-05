%%
%% This is the result of comparison with Zhang's method and Liu's method
%%

%% Definition of some initial variables
clc;
clear;
global NStep NUMC PECV_gammaTrue startGamma;
global Freq Width Height midWidth midHeight sizeStep base Beta Lambda centerROI;
accuracy = 3;

Freq = 64;                      % Frequency of the projected sine diagram during the test
NStep = 3;                      % Phase shift steps of the projected sine diagram during the test
Width = 1920;                   % The width of the projected sine diagram during the test
Height = 1080;                  % The height of the projected sine diagram during the test
Beta = 0.5;                     % Other variables of the projected sine diagram during the test
Lambda = 1 ;                    % Other variables of the projected sine diagram during the test
midHeight = round(Height/2);    % the midpoint of the projected sine diagram on the height during the test
midWidth = round(Width/2);      % the midpoint of the projected sine diagram on the width during the test
sizeStep = 150;                 % A [(2*sizeStep) * (2*sizeStep)]  rectangular ROI is used to solve for the ideal phase
NUMC = 3;                       % Dimensionality of PECV
startGamma = 1;                 % Initial gamma value of the PECV method iteration
centerROI = 100;                % A [centerROI * centerROI]  rectangular ROI is used to solve for the ideal phase

%% find current path
currentPath = mfilename('fullpath');
i=findstr(currentPath,'\');
currentPath=currentPath(1:i(end));
path = [currentPath,'guesspattens'];

savePath = currentPath;
FontName = 'Times New Roman';
FontSize = 12;
FontWeight = 'bold';
saveType = '-dpng';
saveResolution = '-r300'; 


%% Experiment

%% calculate the gamma value 

% zhang
NStep_Zhang = 16;
[B1_f1, B2_f1, B3_f1, filter] = getBifromPictureSetZhang(path, NStep_Zhang, 0);
[B1_f2, B2_f2, B3_f2, filter] = getBifromPictureSetZhang(path, NStep_Zhang, NStep_Zhang);
[B1_f3, B2_f3, B3_f3, filter] = getBifromPictureSetZhang(path, NStep_Zhang, 2*NStep_Zhang);

B2_B1 = B2_f1 ./ B1_f2;
B3_B1 = B3_f1 ./ B1_f3;

B2_B1 = B2_B1.*filter;
B3_B1 = B3_B1.*filter;

% k = 1 because we have chosen B1 B2 B3 to calculate the gamma value
k = 1;
beta = 0.5;
gamma = ((k + 2) * B3_B1 + k - (2*beta - 2)/beta * (k + 1) * B2_B1) ./ (1 - B3_B1);

gamma = gamma.*filter;
gamma(find(isnan(gamma)==1)) = 0;

ZhangGammaValue = mean(gamma(gamma~=0));

% % Liu's algorithm
NStep_Liu = NStep_Zhang;
[B1, B2, filter] = getBifromPictureSetLiu(path, NStep_Liu, 2*NStep_Zhang);
gamma = (B1 + 2 * B2) ./ (B1 - B2);
gamma = gamma.*filter;
gamma(find(isnan(gamma)==1)) = 0;

LiuGammaValue = mean(gamma(gamma~=0));

% our algorithm
NStep_Our = 3;
[PhaseError1, PhaseError2] = getPhase(path, NStep_Our, 3*NStep_Zhang);
PhaseError1 = PhaseError1(1650:1950 ,1800:2100 );
PhaseError2 = PhaseError2(1650:1950 ,1800:2100 );
PhaseError1 = unwrap_phase(PhaseError1);
PhaseError2 = unwrap_phase(PhaseError2);

IdealPhase = getIdealPhasebyTwoPhase(PhaseError1, PhaseError2);
PECV_gammaTrue = getPECV(getCenterROI(PhaseError1),getCenterROI(IdealPhase),NUMC,NStep_Our);
[bestGamma,loss,exitflag,output] = fminsearch(@myPaperfun,startGamma);
PECV_gammaTrue = getPECV(getCenterROI(PhaseError2),getCenterROI(IdealPhase),NUMC,NStep_Our);
[bestGamma2,loss,exitflag,output] = fminsearch(@myPaperfun,startGamma);

gammaA = 2 * (bestGamma - bestGamma2);
gammaB = bestGamma - gammaA;
OurGammaValue = gammaA / (1 - gammaB);


% ZhangGammaValue=2.3489   LiuGammaValue=2.0951  OurGammaValue=2.3238
disp(['ZhangGammaValue=',num2str(ZhangGammaValue), '   LiuGammaValue=', num2str(LiuGammaValue), '  OurGammaValue=', num2str(OurGammaValue)])

%% the results after correction 

path = [currentPath,'finalresults'];
[ZhangPhase, LiuPhase, OurPhase, NoCorrectPhase, IdealPhase] = getAllPhase(path);

% The data of 1000 * 1000 center region of phase were taken for analysis
startWidth = 1100;
endWidth = 3100;
startHeight = 1100;
endHeight = 2100;

ZhangPhase = ZhangPhase(startHeight:endHeight,startWidth:endWidth);
LiuPhase = LiuPhase(startHeight:endHeight,startWidth:endWidth);
OurPhase = OurPhase(startHeight:endHeight,startWidth:endWidth);
NoCorrectPhase = NoCorrectPhase(startHeight:endHeight,startWidth:endWidth);
IdealPhase = IdealPhase(startHeight:endHeight,startWidth:endWidth);

ZhangPhase = unwrap_phase(ZhangPhase);
LiuPhase = unwrap_phase(LiuPhase);
OurPhase = unwrap_phase(OurPhase);
NoCorrectPhase = unwrap_phase(NoCorrectPhase);
IdealPhase = unwrap_phase(IdealPhase);

ERROR_NoCorrect = NoCorrectPhase - IdealPhase;
ERROR_ZhangCorrect = ZhangPhase - IdealPhase;
ERROR_LiuCorrect = LiuPhase - IdealPhase;
ERROR_OURCorrect = OurPhase - IdealPhase;

MaxNoCorrect = max(max(abs(ERROR_NoCorrect)));
MaxZhangCorrect = max(max(abs(ERROR_ZhangCorrect)));
MaxLiuCorrect = max(max(abs(ERROR_LiuCorrect)));
MaxOURCorrect = max(max(abs(ERROR_OURCorrect)));

STDNoCorrect = CalculateSTDERROR(ERROR_NoCorrect);
STDZhangCorrect = CalculateSTDERROR(ERROR_ZhangCorrect);
STDLiuCorrect = CalculateSTDERROR(ERROR_LiuCorrect);
STDOURCorrect = CalculateSTDERROR(ERROR_OURCorrect);

MAENoCorrect = CalculateMAEERROR(ERROR_NoCorrect);
MAEZhangCorrect = CalculateMAEERROR(ERROR_ZhangCorrect);
MAELiuCorrect = CalculateMAEERROR(ERROR_LiuCorrect);
MAEOURCorrect = CalculateMAEERROR(ERROR_OURCorrect);

row = 500;
FF = figure;
plot(ERROR_NoCorrect(row,1:endWidth-startWidth),'g');
hold on
plot(ERROR_ZhangCorrect(row,1:endWidth-startWidth),'c');
hold on
plot(ERROR_LiuCorrect(row,1:endWidth-startWidth),'b');
hold on
plot(ERROR_OURCorrect(row,1:endWidth-startWidth),'r');
hold on
plot(1:endWidth-startWidth,0 * (1:endWidth-startWidth),':','Color','k','linewidth', 1.5);

h=legend({['Without Correction',newline,'STD: ',num2str(STDNoCorrect),newline,'MAX: ',num2str(MaxNoCorrect),newline,'MAE: ',num2str(MAENoCorrect)],
    ['Zhang Method Correction',newline,'STD: ',num2str(STDZhangCorrect),newline,'MAX: ',num2str(MaxZhangCorrect),newline,'MAE: ',num2str(MAEZhangCorrect)],
    ['Liu Method Correction',newline,'STD: ',num2str(STDLiuCorrect),newline,'MAX: ',num2str(MaxLiuCorrect),newline,'MAE: ',num2str(MAELiuCorrect)],
    ['Our Method Correction',newline,'STD: ',num2str(STDOURCorrect),newline,'MAX: ',num2str(MaxOURCorrect),newline,'MAE: ',num2str(MAEOURCorrect)]
    },'Location','north','NumColumns',4);
set(h,'FontName',FontName,'FontSize',FontSize-3,'FontWeight',FontWeight) 
ylabel(' Phase Error(rad) ');
xlabel(' Pixel ');
set(gca,'FontName',FontName,'FontSize',FontSize, 'FontWeight',FontWeight);
set(gca,'YLim',[-MaxNoCorrect-0.01 MaxNoCorrect+0.2]);
set(gca,'XTickLabel',[startWidth:200:endWidth]);
w = get(gca,'position'); 
original = get(FF,'position');
set(FF,'position',[original(1)*0.5, original(2)*0.5, original(3)*2.5, original(4)*1.3]);
% print(gcf,[savePath,'PhaseError'],saveType,saveResolution);
% close

%%  frequency spectrum
path = [currentPath,'finalresults'];

ZhangCorrect = double( imread( sprintf( '%s/%02d.bmp', path, 0 )));
LiuCorrect = double( imread( sprintf( '%s/%02d.bmp', path, 3 )));
OurCorrect = double( imread( sprintf( '%s/%02d.bmp', path, 6 )));
NoCorrect = double( imread( sprintf( '%s/%02d.bmp', path, 9 )));

subplot(2,2,1);
f1=fftshift(abs(fft(NoCorrect(2000,:))),2); 
plot(f1(1700:2300),'r-')

    ylabel(' Amplitude ');
    xlabel(' Frequency ');
    set(gca,'XTickLabel',[1800:50:2200]);

subplot(2,2,2);
f2=fftshift(abs(fft(ZhangCorrect(2000,:))),2);
plot(f2(1700:2300),'r-')

    ylabel(' Amplitude ');
    xlabel(' Frequency ');
    set(gca,'XTickLabel',[1800:50:2200]);


subplot(2,2,3);
f3=fftshift(abs(fft(LiuCorrect(2000,:))),2);
plot(f3(1700:2300),'r-')

    ylabel(' Amplitude ');
    xlabel(' Frequency ');
    set(gca,'XTickLabel',[1800:50:2200]);

subplot(2,2,4);
f4=fftshift(abs(fft(OurCorrect(2000,:))),2);
plot(f4(1700:2300),'r-')

    ylabel(' Amplitude ');
    xlabel(' Frequency ');
    set(gca,'XTickLabel',[1800:150:2200]);

close;

% f = {f1,f2,f3,f4};
% saveName = {'NoCorrect','ZhangCorrect','LiuCorrect','OurCorrect'};
% % save
% for i = 1 : length(f)
%     plot(f{i},'r-');
% %     set(gca,'XLim',[-50 2000]);
%     ylabel(' Amplitude ');
%     xlabel(' Frequency ');
%     print(gcf,[savePath,saveName{i}],saveType,saveResolution);
%     close
% end
% 
% % Detail View
% for i = 1 : length(f)
%     plot(f{i}(1700:2300),'r-');
%     ylabel(' Amplitude ');
%     xlabel(' Frequency ');
%     set(gca,'YLim',[-10 20000]);
%     set(gca,'XTickLabel',[1800:50:2200]);
%     print(gcf,[savePath,saveName{i},'big'],saveType,saveResolution);
%     close
% end



%% Simulation
Freq_Zhang = 22;
Freq_Liu = 64;
Freq_Our = 64;
NStep_Zhang = 16;
NStep_Liu = 16;
NStep_Our = 3;

global PictureSet_Zhang PictureSet_Liu;
PictureSet_Zhang = cell(3, NStep_Zhang);
PictureSet_Liu = cell(1, NStep_Liu);

for i = 1 : NStep_Zhang
    n = i-1;
    Im1 = Picture(Freq_Zhang, NStep_Zhang,1920,1080,n,1,0.5,1);
    Im2 = Picture(2 * Freq_Zhang, NStep_Zhang,1920,1080,n,1,0.5,1);
    Im3 = Picture(3 * Freq_Zhang, NStep_Zhang,1920,1080,n,1,0.5,1);
    PictureSet_Zhang{1,i} = Im1;
    PictureSet_Zhang{2,i} = Im2;
    PictureSet_Zhang{3,i} = Im3;
end
  
for i = 1 : NStep_Liu
    n = i-1;
    Im = Picture(Freq_Liu, NStep_Liu,1920,1080,n,1,0.5,1);
    PictureSet_Liu{1,i} = Im;
end


% save the data
STDGammaError = zeros(10,4);
STDGammaError_midcase = zeros(201,4);
Index = 1;
Index2 = 1; 

% impactFactor range
GaussianNoise = 0:0.0001:0.001; 

% choose impactFactor
impactFactor = 'Noise'; impactFactorRange = GaussianNoise;

for factor = impactFactorRange 
ZhangLossAll = 0;
LiuLossAll = 0;
OurLossAll = 0;
    
    for testGamma = 1:0.01:3
        
        % Zhang's algorithm
        startTime = clock; 
        [B1_f1, B2_f1, B3_f1, filter] = getBifromPictureSetZhang_simulation(1, NStep_Zhang, testGamma, factor, impactFactor);
        [B1_f2, B2_f2, B3_f2, filter] = getBifromPictureSetZhang_simulation(2, NStep_Zhang, testGamma, factor, impactFactor);
        [B1_f3, B2_f3, B3_f3, filter] = getBifromPictureSetZhang_simulation(3, NStep_Zhang, testGamma, factor, impactFactor);

        B2_B1 = B2_f1 ./ B1_f2;
        B3_B1 = B3_f1 ./ B1_f3;

        B2_B1 = B2_B1;
        B3_B1 = B3_B1;

        % k = 1 because we have chosen B1 B2 B3 to calculate the gamma value
        k = 1;
        beta = 0.5;
        gamma = ((k + 2) * B3_B1 + k - (2*beta - 2)/beta * (k + 1) * B2_B1) ./ (1 - B3_B1);

        ZhangGammaValue = mean(gamma(gamma~=0));
        ZhangGammaError = abs(ZhangGammaValue - testGamma);
        ZhangLossAll = ZhangLossAll + ZhangGammaError.^2;
        if factor == impactFactorRange(round(length(impactFactorRange)/2,0))
            STDGammaError_midcase(Index2,1) = testGamma;
            STDGammaError_midcase(Index2,2) = ZhangGammaError;
        end

        endTime = clock;
        
    disp([' factor / gamma  ', num2str(factor),' / ', num2str(testGamma)  ,'          ','ZhangGammaValue=',num2str(ZhangGammaValue), '    ZhangGammaError=', num2str(ZhangGammaError) , '    time-cost=' , num2str(etime(endTime,startTime))]);

%         Liu's algorithm
        startTime = clock;
        [B1, B2, filter] = getBifromPictureSetLiu_simulation(Freq_Liu, NStep_Liu, testGamma, factor, impactFactor);
        gamma = (B1 + 2 * B2) ./ (B1 - B2);
        LiuGammaValue = mean(gamma(gamma~=0));
        LiuGammaError = abs(LiuGammaValue - testGamma);
        LiuLossAll = LiuLossAll + LiuGammaError.^2;
        if factor == impactFactorRange(round(length(impactFactorRange)/2,0))
            STDGammaError_midcase(Index2,3) = LiuGammaError;
        end
        endTime = clock;
     disp([' factor / gamma  ', num2str(factor),' / ', num2str(testGamma)  ,'          ','LiuGammaValue=',num2str(LiuGammaValue), '    LiuGammaError=', num2str(LiuGammaError) , '    time-cost=' , num2str(etime(endTime,startTime))]);


        % our algorithm
        startTime=clock;

        [PhaseError1, PhaseError2] = getPhase_simulation(Freq_Our, NStep_Our, testGamma, factor, impactFactor);
        PhaseError1 = unwrap_phase(PhaseError1); % PhaseUnwrapping
        PhaseError2 = unwrap_phase(PhaseError2);
        IdealPhase = getIdealPhasebyTwoPhase(PhaseError1, PhaseError2);
        
       startTime=clock;
        PECV_gammaTrue = getPECV(PhaseError1,IdealPhase,NUMC,NStep);
        [bestGamma,loss,exitflag,output] = fminsearch(@myPaperfun,startGamma); 
        PECV_gammaTrue = getPECV(PhaseError2,IdealPhase,NUMC,NStep); % getCenterROI
        [bestGamma2,loss,exitflag,output] = fminsearch(@myPaperfun,startGamma);

        gammaA = 2 * (bestGamma - bestGamma2);
        gammaB = bestGamma - gammaA;
        OurGammaValue = gammaA / (1 - gammaB);
        
        endTime = clock;
        
        OurGammaError = abs(OurGammaValue - testGamma);
        OurLossAll = OurLossAll + OurGammaError.^2;

        if factor == impactFactorRange(round(length(impactFactorRange)/2,0))
            STDGammaError_midcase(Index2,4) = OurGammaError;
            Index2 = Index2 + 1;
        end

    disp([' factor / gamma  ', num2str(factor),' / ', num2str(testGamma)  ,'          ','OurGammaValue=',num2str(OurGammaValue), '    OurGammaError=', num2str(OurGammaError)  ,'    time-cost=' , num2str(etime(endTime,startTime))]);

    end

    STDGammaError(Index,1) = factor;
    STDGammaError(Index,2) = sqrt(ZhangLossAll)/201; % Here are 201 tests from 1 to 3 increasing by 0.01 each time
    STDGammaError(Index,3) = sqrt(LiuLossAll)/201;
    STDGammaError(Index,4) = sqrt(OurLossAll)/201;
    Index = Index + 1;

end

save STDGammaErrorZLO.mat  STDGammaError;
save STDGammaError_midcase.mat  STDGammaError_midcase;


load('STDGammaErrorZLO.mat')

noise = STDGammaError(:,1);
ZhangMethod = STDGammaError(:,2);
LiuMethod = STDGammaError(:,3);
OurMethod = STDGammaError(:,4);

plot(noise, ZhangMethod,'r', 'linewidth', 1.5);
hold on
plot(noise, LiuMethod,'g', 'linewidth', 1.5);
hold on
plot(noise, OurMethod,'b', 'linewidth', 1.5);

h=legend({'Zhang Method','Liu Method','Our Method'},'Location','northwest','NumColumns',1);
set(h,'FontName',FontName,'FontSize',FontSize-1,'FontWeight',FontWeight)
ylabel(' STDGammaError');
xlabel(' Gaussian Noise Level');
set(gca,'FontName',FontName,'FontSize',FontSize, 'FontWeight',FontWeight);
% grid on;
print(gcf,[savePath,'NoiseSTD'],saveType,saveResolution)
close


% 0.0005 noise
load('STDGammaError_midcase.mat')

Gamma = STDGammaError_midcase(:,1);
ZhangMethod = STDGammaError_midcase(:,2);
LiuMethod = STDGammaError_midcase(:,3);
OurMethod = STDGammaError_midcase(:,4);

plot(Gamma, ZhangMethod,'r', 'linewidth', 1.5);
hold on
plot(Gamma, LiuMethod,'g', 'linewidth', 1.5);
hold on
plot(Gamma, OurMethod,'b', 'linewidth', 1.5);
h=legend({'Zhang Method','Liu Method','Our Method'},'Location','northwest','NumColumns',1);
set(h,'FontName',FontName,'FontSize',FontSize-1,'FontWeight',FontWeight)
ylabel(' absGammaError ');
xlabel(' Gamma ');
set(gca,'FontName',FontName,'FontSize',FontSize, 'FontWeight',FontWeight);
set(gca,'YLim',[-0.005 0.35]);
% grid on;
print(gcf,[savePath,'0005NoiseSTD'],saveType,saveResolution)
close



%% This is a test of the monotonically increasing nature of the sum of the elements of PECV that we exploit (in the case beta = 0.5 lambda = 1)


Beta = 0.5;
Lambda = 1 ;
base = createPhaseImage(1,NStep,Width,1,1,Beta,Lambda);

% To test PECV in higher dimensions, width must be increased
for N = 3 : 150
    SUM_OF_PECV_NUMC = [];
    
    INDEX2 = 1;
    for gamma = 1:0.01:3
        ErrorPhase = createPhaseImage(1,NStep,Width,1,gamma,Beta,Lambda);
        PECV_gammaTrue = getPECV(ErrorPhase,base,N,NStep);
        SUM_OF_PECV_NUMC = [SUM_OF_PECV_NUMC, sum(PECV_gammaTrue)];
    end
    
    if IsMonotoneIncreasing(SUM_OF_PECV_NUMC) == 1
        str = [ 'N = ',num2str(N) , ' IsMonotoneIncreasing:  [True]'];
    else
        str = [ 'N = ',num2str(N) , ' IsMonotoneIncreasing:  [False]'];
    end
    
    disp(str)
end






%% Zhang's model
function [B1, B2, B3, filter] = getBifromPictureSetZhang(path, NStep, Count)
    I1 = 0; I2 = 0; I3 = 0;
    PSin = 0;   PCos = 0;

    for i = 1 : NStep
        n = i-1;
        Im = double( imread( sprintf( '%s/%02d.bmp', path, Count )));
        Count = Count + 1;
        
        I1 = I1 + Im*exp(-1i * 2*pi * n * 1 / NStep);
        I2 = I2 + Im*exp(-1i * 2*pi * n * 2 / NStep);
        I3 = I3 + Im*exp(-1i * 2*pi * n * 3 / NStep);
        
        PSin = PSin + Im * sin(2*pi * n / NStep);
        PCos = PCos + Im * cos(2*pi * n / NStep);
    end

    Bc = 2 / NStep * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    filter = Index;

    B1 = 2/NStep * abs(I1);
    B2 = 2/NStep * abs(I2);
    B3 = 2/NStep * abs(I3);

    if( abs(angle(I2) - 2*angle(I1)) > pi/2 )
        B2 = -B2;
    end

    if( abs(angle(I3) - 3*angle(I1)) > pi/2 )
        B3 = -B3;
    end
end

%% Liu's model
function [B1, B2, filter] = getBifromPictureSetLiu(path, NStep, Count)
    PSin1 = 0;   PCos1 = 0;
    PSin2 = 0;   PCos2 = 0;
    
    for i = 1 : NStep
        n = i-1;
        Im = double( imread( sprintf( '%s/%02d.bmp', path, Count )));
        Count = Count + 1;
       
        PSin1 = PSin1 + Im * sin(1 * 2*pi * n / NStep);
        PCos1 = PCos1 + Im * cos(1 * 2*pi * n / NStep);
        PSin2 = PSin2 + Im * sin(2 * 2*pi * n / NStep);
        PCos2 = PCos2 + Im * cos(2 * 2*pi * n / NStep);
    end

    Bc = 2 / NStep * sqrt( PSin1 .* PSin1 + PCos1 .* PCos1 );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    filter = Index;

    B1 = 2/NStep * sqrt( PSin1.*PSin1 +  PCos1 .* PCos1 );
    B2 = 2/NStep * sqrt( PSin2.*PSin2 +  PCos2 .* PCos2 );
end


function [PhaseError1, PhaseError2] = getPhase(path, NStep,Count)

    PSin1 = 0;   PCos1 = 0;
    PSin2 = 0;   PCos2 = 0;
    global W;
    
    for i = 1 : NStep
        n = i-1;
        Im = double( imread( sprintf( '%s/%02d.bmp', path, Count )));
        Im2 = double( imread( sprintf( '%s/%02d.bmp', path, Count+NStep )));
        Count = Count + 1;
    
        PSin1 = PSin1 + Im * sin(2*pi * n / NStep);
        PCos1 = PCos1 + Im * cos(2*pi * n / NStep);
        PSin2 = PSin2 + Im2 * sin(2*pi * n / NStep);
        PCos2 = PCos2 + Im2 * cos(2*pi * n / NStep);
    end

    Bc = 2 / NStep * sqrt( PSin1 .* PSin1 + PCos1 .* PCos1 );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;

    PhaseError1 = pi + atan2(PSin1,  -PCos1 );
    PhaseError2 = pi + atan2(PSin2,  -PCos2 );
    PhaseError1 = PhaseError1.*Index;
    PhaseError2 = PhaseError2.*Index;
end

function [ZhangPhase, LiuPhase, OurPhase, NoCorrectPhase, IdealPhase] = getAllPhase(path)
    
    I0 = double(imread( sprintf( '%s/%02d.bmp', path, 0 )) );
    I1 = double(imread( sprintf( '%s/%02d.bmp', path, 1 )) );
    I2 = double(imread( sprintf( '%s/%02d.bmp', path, 2 )) );
 
    PSin = I0 * sin(2*pi * 0 / 3) + I1 * sin(2*pi * 1 / 3) + I2 * sin(2*pi * 2 / 3);
    PCos = I0 * cos(2*pi * 0 / 3) + I1 * cos(2*pi * 1 / 3) + I2 * cos(2*pi * 2 / 3);

    Bc = 2 / 3 * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    
    ZhangPhase = (pi + atan2(  PSin,-PCos )).*Index;

    
    I0 = double(imread( sprintf( '%s/%02d.bmp', path, 3 )) );
    I1 = double(imread( sprintf( '%s/%02d.bmp', path, 4 )) );
    I2 = double(imread( sprintf( '%s/%02d.bmp', path, 5 )) );
 
    PSin = I0 * sin(2*pi * 0 / 3) + I1 * sin(2*pi * 1 / 3) + I2 * sin(2*pi * 2 / 3);
    PCos = I0 * cos(2*pi * 0 / 3) + I1 * cos(2*pi * 1 / 3) + I2 * cos(2*pi * 2 / 3);

    Bc = 2 / 3 * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    
    LiuPhase = (pi + atan2(  PSin,-PCos )).*Index;
    
    I0 = double(imread( sprintf( '%s/%02d.bmp', path, 6 )) );
    I1 = double(imread( sprintf( '%s/%02d.bmp', path, 7 )) );
    I2 = double(imread( sprintf( '%s/%02d.bmp', path, 8 )) );
 
    PSin = I0 * sin(2*pi * 0 / 3) + I1 * sin(2*pi * 1 / 3) + I2 * sin(2*pi * 2 / 3);
    PCos = I0 * cos(2*pi * 0 / 3) + I1 * cos(2*pi * 1 / 3) + I2 * cos(2*pi * 2 / 3);

    Bc = 2 / 3 * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    
    OurPhase = (pi + atan2(  PSin,-PCos )).*Index;
    
    
    I0 = double(imread( sprintf( '%s/%02d.bmp', path, 9 )) );
    I1 = double(imread( sprintf( '%s/%02d.bmp', path, 9 +12 )) );
    I2 = double(imread( sprintf( '%s/%02d.bmp', path, 9 + 24 )) );
 
    PSin = I0 * sin(2*pi * 0 / 3) + I1 * sin(2*pi * 1 / 3) + I2 * sin(2*pi * 2 / 3);
    PCos = I0 * cos(2*pi * 0 / 3) + I1 * cos(2*pi * 1 / 3) + I2 * cos(2*pi * 2 / 3);

    Bc = 2 / 3 * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    
    NoCorrectPhase = (pi + atan2(  PSin,-PCos )).*Index;
    
    
    PSin = 0;   PCos = 0; NStep = 36;
    Count = 9;
    for i = 1 : NStep
        n = i-1;
        Im = double( imread( sprintf( '%s/%02d.bmp', path, Count )));
        Count = Count + 1;

        PSin = PSin + Im * sin(2*pi * n / NStep);
        PCos = PCos + Im * cos(2*pi * n / NStep);      
    end

    Bc = 2 / NStep * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    
    IdealPhase = pi + atan2(PSin,  -PCos ).*Index;
end













%% PECV method
% getPECV
function [PECV] = getPECV(phaseError,phaseTure,NUMC,N)
    [height,width] = size(phaseError);
    NUM = width * height;

    A = zeros(NUM,NUMC);
    B = zeros(NUM,1);
    for m = 1 : NUM
        row = ceil(m/width);
        col = m - (row-1) * width;

        for n = 1 : NUMC  
            A(m,n) = sin(n*N * phaseTure(row,col));
        end

        B(m,1) = phaseTure(row,col) - phaseError(row,col);
    end

    PECV = (A'*A)\(A'*B);
end

% calculate the PECV corresponding to the given gamma value in the case of N-step phase shift 
function [PECV] = CalculatePECV_ByGamma(Gamma,NUMC,testNstep)
    global Beta Lambda;
    width = 100;
    height = 1;

    phaseError = createPhaseImage(1,testNstep,width,height,Gamma, Beta, Lambda);
    phaseTure = createPhaseImage(1,testNstep,width,height,1, Beta, Lambda);
    
    PECV = getPECV(phaseError, phaseTure, NUMC, testNstep);
end

function [centerPicture] = getCenterROI(phase)
    global centerROI;
    [height,width] = size(phase);
    startWidth = round((width - centerROI)/2);
    startHeight = round((height - centerROI)/2);
    centerPicture = phase(startHeight : startHeight + centerROI, startWidth : startWidth + centerROI);
end

function similarity = myPaperfun(gamma)
    global PECV_gammaTrue;
    global NStep NUMC;
    PECV_gammaGuess = CalculatePECV_ByGamma(gamma,NUMC,NStep);
    temp = (PECV_gammaGuess - PECV_gammaTrue).^2;
    similarity = mean(temp(:));
end

%% get ideal phase

function [IdealPhase] = getIdealPhasebyTwoPhase(phaseError1, phaseError2)
    global NUMC NStep;
    NUMC = 3;
    
    x = 1:size(phaseError1,2);
    t=polyfit(x,phaseError1(1,:),1);
    y = t(1)*x+t(2);
    a = sum(abs(phaseError1(1,:) - y));
    b = sum(abs(phaseError2(1,:) - y));
    
    if( a <= b)
        phaseError = phaseError1;
    else
        phaseError = phaseError2;
    end

    [height,width] = size(phaseError);
    phaseTure = findPlane(phaseError);
    
    

    PECV_first = getPECV(getCenterROI2(phaseError), getCenterROI2(phaseTure), NUMC, NStep);
    phaseError2 = CompensateC(phaseError,NStep, PECV_first); 

    IdealPhase = findPlane(phaseError2);
end

function [P] = getCenterROI2(Image)
    [height,width] = size(Image);
    
    length1 = round(width/10);
    start = round((width - length1)/2);
    
    P = Image(:, start:start + length1);
    
    
end



% using the PECV compensate the error phase
function [BestPhase] = CompensateC(Image, NStep, vectorC)
    numC = length(vectorC);

    step = 0;
    Xk = Image;
    while 1
        Xk1 = Image;
        for i = 1:numC
            Xk1 = Xk1 + vectorC(i) * sin(i*NStep*Xk);
        end
        temp = Xk1 - Xk;
        step = step + 1;
        if abs(max(temp(:)))<0.0001 || step >= 30 %收敛判断
            break;
        end
        Xk = Xk1;
    end
    BestPhase = Xk;
end

function [P] = findPlane(Image)
    width = size(Image,2);
    height = size(Image,1);
    N = width*height;
 
    X = height*((1+width)*width/2);
    Y = width*((1+height)*height/2);
    XY = ((1+width)*width/2)*((1+height)*height/2);
    X_2 = height*(width*(width+1)*(2*width+1)/6);
    Y_2 = width*(height*(height+1)*(2*height+1)/6);
    
    A = zeros(3,3);
    A(1,1) = X_2;A(2,2) = Y_2;A(3,3) = N;
    A(1,2) = XY;A(2,1) = XY; 
    A(1,3) = X;A(3,1) = X;
    A(2,3) = Y;A(3,2) = Y;
    
    B = zeros(3,1);
    for yi = 1:height
        for xi = 1:width
            B(1,1) = B(1,1) + xi *Image(yi,xi);
            B(2,1) = B(2,1) + yi *Image(yi,xi);
            B(3,1) = B(3,1) + Image(yi,xi);
        end
    end
    
    C = lsqminnorm(A,B);
    
    P = zeros(height,width);
    for i = 1:height
        for j = 1:width 
            P(i,j) = C(1,1) * j + C(2,1) * i + C(3,1);
        end
    end
end



%% Auxiliary Functions
function [PhaseImage] = createPhaseImage(Freq,Steps,Width,Height,randomGamma,beta,lambda)
    PSin = 0; PCos = 0;
    for step = 1:Steps
        n = step - 1;
        
        Im = Picture(Freq,Steps,Width,Height,n,randomGamma,beta,lambda); 
        PSin = PSin + Im * sin(  2 * pi * n / Steps );
        PCos = PCos + Im * cos(  2 * pi * n / Steps );
        clear Im;
    end
    
    PhaseImage = pi + atan2(PSin,  -PCos );
    PhaseImage(:,1) = 0; 
end

function [P] = Picture(freq,Nset,width,height,index,gamma,beta,lambda)
    Im = zeros(height,width);
    for i = 1:width 
        Im(:,i) = beta + (1 - beta) * lambda * cos(2*pi*(i-1) * freq /width + index*2*pi/Nset);
    end
    P = Im .^(gamma);
end

function [IsTrue] = IsMonotoneIncreasing(TMP)
    Now = TMP(1);
    IsTrue = 1;
    for index = 2:length(TMP)
        if(Now >= TMP(index))
            IsTrue = 0;
            break;
        end
        Now = TMP(index);
    end
end

function [P] = PhaseUnwrapping(Image)
    global Freq NStep Width Height base;
    %base = createPhaseImage(1,NStep,Width,Height,1);
    K = round((Freq*base - Image)/(2*pi));
    P = Image + 2 * K * pi;
end


function [STDERROR] = CalculateSTDERROR(Error)
    [height,width] = size(Error);
    Count = 0;
    ALLERROR = 0;
    for y = 1:height
        for x = 1:width-1
            if(Error(y,x) == 0)  continue; end
            
            Count = Count + 1;
            ALLERROR = ALLERROR + (Error(y,x)-0)^2;
        end
    end
    STDERROR = sqrt(ALLERROR)/sqrt(Count);
end

function [STDERROR] = CalculateMAEERROR(Error)
    [height,width] = size(Error);
    Count = 0;
    ALLERROR = 0;
    for y = 1:height
        for x = 1:width-1
            if(Error(y,x) == 0)  continue; end
            
            Count = Count + 1;
            ALLERROR = ALLERROR + abs(Error(y,x)-0);
        end
    end
    STDERROR = (ALLERROR)/(Count);
end

%% The simulation function

% zhang's model
function [B1, B2, B3, filter] = getBifromPictureSetZhang_simulation(Freq, NStep, gamma, var, select)
    I1 = 0; I2 = 0; I3 = 0;
    PSin = 0;   PCos = 0;
    global PictureSet_Zhang;
    
    for i = 1 : NStep
        n = i-1;
%         Im = Picture(Freq, NStep,1920,1080,n,gamma,0.5,1);
        Im = PictureSet_Zhang{Freq,i}.^gamma;
        
        if strcmp(select,'Noise')
            Im = imnoise(Im,'gaussian',0,var);
        elseif strcmp(select,'Defocusing')
            Im = imfilter(Im, W, 'replicate');
        else
            Im = Im;
        end

        I1 = I1 + Im*exp(-1i * 2*pi * n * 1 / NStep);
        I2 = I2 + Im*exp(-1i * 2*pi * n * 2 / NStep);
        I3 = I3 + Im*exp(-1i * 2*pi * n * 3 / NStep);
        
        PSin = PSin + Im * sin(2*pi * n / NStep);
        PCos = PCos + Im * cos(2*pi * n / NStep);
    end

    Bc = 2 / NStep * sqrt( PSin .* PSin + PCos .* PCos );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    filter = Index;

    B1 = 2/NStep * abs(I1);
    B2 = 2/NStep * abs(I2);
    B3 = 2/NStep * abs(I3);

    if( abs(angle(I2) - 2*angle(I1)) > pi/2 )
        B2 = -B2;
    end

    if( abs(angle(I3) - 3*angle(I1)) > pi/2 )
        B3 = -B3;
    end

end

% Liu's model
function [B1, B2, filter] = getBifromPictureSetLiu_simulation(Freq, NStep, gamma, var, select)
    PSin1 = 0;   PCos1 = 0;
    PSin2 = 0;   PCos2 = 0;
    global PictureSet_Liu;
    
    for i = 1 : NStep
        n = i-1;
%         Im = Picture(Freq, NStep,1920,1080,n,gamma,0.5,1); 
        Im = PictureSet_Liu{1,i}.^gamma;

        if strcmp(select,'Noise')
            Im = imnoise(Im,'gaussian',0,var);
        end
        
        PSin1 = PSin1 + Im * sin(1 * 2*pi * n / NStep);
        PCos1 = PCos1 + Im * cos(1 * 2*pi * n / NStep);
        PSin2 = PSin2 + Im * sin(2 * 2*pi * n / NStep);
        PCos2 = PCos2 + Im * cos(2 * 2*pi * n / NStep);
    end

    Bc = 2 / NStep * sqrt( PSin1 .* PSin1 + PCos1 .* PCos1 );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    filter = Index;

    B1 = 2/NStep * sqrt( PSin1.*PSin1 +  PCos1 .* PCos1 );
    B2 = 2/NStep * sqrt( PSin2.*PSin2 +  PCos2 .* PCos2 );

end

% our model
% Here it is necessary that PhaseError1 has no pre-set gamma value and PhaseError2 has a pre-set gamma value
function [PhaseError1, PhaseError2] = getPhase_simulation(Freq, NStep, gamma, var, select)
    PSin1 = 0;   PCos1 = 0;
    PSin2 = 0;   PCos2 = 0;
    
    for i = 1 : NStep
        n = i-1;
        Im = Picture(Freq,NStep,90000, 1,n,gamma,0.5,1);
        Im2 = Picture(Freq,NStep,90000, 1,n,gamma/2,0.5,1);
        if strcmp(select,'Noise')
            Im = imnoise(Im,'gaussian',0,var);
            Im2 = imnoise(Im2,'gaussian',0,var);
        end
     
        PSin1 = PSin1 + Im * sin(2*pi * n / NStep);
        PCos1 = PCos1 + Im * cos(2*pi * n / NStep);
        PSin2 = PSin2 + Im2 * sin(2*pi * n / NStep);
        PCos2 = PCos2 + Im2 * cos(2*pi * n / NStep);
    end

    Bc = 2 / NStep * sqrt( PSin1 .* PSin1 + PCos1 .* PCos1 );
    Index = zeros( size( Bc ) );
    Idx = find( Bc > 10 );
    Index( Idx ) = 1;
    
    PhaseError1 = pi + atan2(PSin1,  -PCos1 );
    PhaseError2 = pi + atan2(PSin2,  -PCos2 );
    PhaseError1(:,1) = 0;
    PhaseError2(:,1) = 0;
end
