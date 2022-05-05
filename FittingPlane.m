clc;
clear;
%%
global NStep NUMC PECV_gammaTrue startGamma;
global Freq Width Height midWidth midHeight sizeStep base Beta Lambda centerROI;
Freq = 64;                      % Frequency of the projected sine diagram during the test
NStep = 3;                      % Phase shift steps of the projected sine diagram during the test
sizeStep = 150;                 % A [(2*sizeStep) * (2*sizeStep)]  rectangular ROI is used to solve for the ideal phase
NUMC = 3;                       % Dimensionality of PECV
startGamma = 1;                 % Initial gamma value of the PECV method iteration
centerROI = 100;
Beta = 0.5;
Lambda = 1;


%% find current path
currentPath = mfilename('fullpath');
i=findstr(currentPath,'\');
currentPath=currentPath(1:i(end));
path = [currentPath,'fitpatterns'];

savePath = currentPath;
FontName = 'Times New Roman';
FontSize = 12;
FontWeight = 'bold';
saveType = '-dpng';
saveResolution = '-r300'; 

%% getting phase
PSin = 0;
PCos = 0;
PSIN1 = 0;
PCOS1 = 0;
Count = 0;
N = 36;
for j = 1 : N
    n = j - 1;           
    Im = double(imread( sprintf( '%s/%02d.bmp', path, Count )) );
    Count = Count + 1;
    PSin = PSin + Im * sin(  2 * pi * n / N );
    PCos = PCos + Im * cos(  2 * pi * n / N );

    if mod(n,2) == 0
        PSIN1 = PSIN1 + Im * sin(  2 * pi * n / N );
        PCOS1 = PCOS1 + Im * cos(  2 * pi * n / N );
    end
    clear Im;
end   

Bc = 2 / N * sqrt( PSin .* PSin + PCos .* PCos );

Index = zeros( size( Bc ) );
Idx = find( Bc > 10 );
Index( Idx ) = 1;
Bc = Bc .* Index;
    
IdealPhase_36step = pi + atan2(  PSin,  -PCos ) .* Index;
IdealPhase_18step = pi + atan2(  PSIN1,  -PCOS1 ) .* Index;

I0 = double(imread( sprintf( '%s/%02d.bmp', path, Count )) );
I1 = double(imread( sprintf( '%s/%02d.bmp', path, Count+1 )) );
I2 = double(imread( sprintf( '%s/%02d.bmp', path, Count+2 )) );
PSIN = I0 * sin(2*pi * 0 / 3) + I1 * sin(2*pi * 1 / 3) + I2 * sin(2*pi * 2 / 3);
PCOS = I0 * cos(2*pi * 0 / 3) + I1 * cos(2*pi * 1 / 3) + I2 * cos(2*pi * 2 / 3);
ErrorPhase2 = (pi + atan2(  PSIN,-PCOS )).*Index;

Count = 0;
I0 = double(imread( sprintf( '%s/%02d.bmp', path, Count )) );
I1 = double(imread( sprintf( '%s/%02d.bmp', path, Count+12 )) );
I2 = double(imread( sprintf( '%s/%02d.bmp', path, Count+24 )) );
PSIN = I0 * sin(2*pi * 0 / 3) + I1 * sin(2*pi * 1 / 3) + I2 * sin(2*pi * 2 / 3);
PCOS = I0 * cos(2*pi * 0 / 3) + I1 * cos(2*pi * 1 / 3) + I2 * cos(2*pi * 2 / 3);
ErrorPhase1 = (pi + atan2(  PSIN,-PCOS )).*Index;

% get center area
startWidth = 1700;
endWidth = 2000;
startHeight = 1700;
endHeight = 2000;
IdealPhase_36step = IdealPhase_36step(startHeight:endHeight,startWidth:endWidth);
IdealPhase_18step = IdealPhase_18step(startHeight:endHeight,startWidth:endWidth);
ErrorPhase2 = ErrorPhase2(startHeight:endHeight,startWidth:endWidth);
ErrorPhase1 = ErrorPhase1(startHeight:endHeight,startWidth:endWidth);

% unwrap phase
IdealPhase_36step = unwrap_phase(IdealPhase_36step);
IdealPhase_18step = unwrap_phase(IdealPhase_18step);
ErrorPhase2 = unwrap_phase(ErrorPhase2);
ErrorPhase1 = unwrap_phase(ErrorPhase1);

%% fitting plane to get ideal phase
IdealPhaseByFitting =  getIdealPhasebyTwoPhase(ErrorPhase1,ErrorPhase2);


%% guess gamma values by three ideal phases
guessGamma_36step = estimateGammaValue(ErrorPhase1, ErrorPhase2, IdealPhase_36step);
guessGamma_18step = estimateGammaValue(ErrorPhase1, ErrorPhase2, IdealPhase_18step);
guessGamma_our = estimateGammaValue(ErrorPhase1, ErrorPhase2, IdealPhaseByFitting);
disp(['guessGamma_36step=', num2str(guessGamma_36step), '     guessGamma_18step=',  num2str(guessGamma_18step),  '     guessGamma_our=', num2str(guessGamma_our)]);

%% residual phase error comparison
ERROR_NoCorrect = getCenterROI(IdealPhase_18step - IdealPhase_36step);           
ERROR_OurCorrect = getCenterROI(IdealPhaseByFitting - IdealPhase_36step);

MaxNoCorrect = max(max(abs(ERROR_NoCorrect)));
MaxOurCorrect = max(max(abs(ERROR_OurCorrect)));

STDNoCorrect = CalculateSTDERROR(ERROR_NoCorrect);
STDOurCorrect = CalculateSTDERROR(ERROR_OurCorrect);

MAENoCorrect = CalculateMAEERROR(ERROR_NoCorrect);
MAEOurCorrect = CalculateMAEERROR(ERROR_OurCorrect);


%% Here, excessive MAX value also indicates excessive noise in the system, because the MAX of 3-steps phase-shifting algorithm is large
%% but the MAX of 18-stepsphase-shifting algorithm is small. 
%% The larger the phase-shifting steps, the higher the anti-noise performance

FF = figure;
plot(ERROR_NoCorrect(78,1:100),'g');
hold on
plot(ERROR_OurCorrect(78,1:100),'r');

h=legend({['18 step',newline,'STD: ',num2str(STDNoCorrect),newline,'MAX: ',num2str(MaxNoCorrect),newline,'MAE: ',num2str(MAENoCorrect)],
    ['Our fitting Method',newline,'STD: ',num2str(STDOurCorrect),newline,'MAX: ',num2str(MaxOurCorrect),newline,'MAE: ',num2str(MAEOurCorrect)]
    },'Location','north','NumColumns',2);
set(h,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight) 
ylabel(' Phase Error(rad) ');
xlabel(' Pixel ');
set(gca,'FontName',FontName,'FontSize',FontSize, 'FontWeight',FontWeight);
set(gca,'YLim',[-0.03 +0.04]);
set(gca,'XTickLabel',[startWidth:30:endWidth]);
w = get(gca,'position'); 
original = get(FF,'position');
set(FF,'position',[original(1)*0.5, original(2)*0.5, original(3)*2.5, original(4)]);
grid on
print(gcf,[savePath,'FitError'],saveType,saveResolution);




%% PECV method
function [finalGamma, timecost, loss, exitflag, output] = estimateGammaValue(phaseError1, phaseError2, IdealPhase)
    global NStep NUMC PECV_gammaTrue startGamma;
    
    % This is the solving part of the PECV algorithm
    startTime=clock;
    PECV_gammaTrue = getPECV(getCenterROI(phaseError1),getCenterROI(IdealPhase),NUMC,NStep);
    [bestGamma,loss,exitflag,output] = fminsearch(@myPaperfun,startGamma);
    PECV_gammaTrue = getPECV(getCenterROI(phaseError2),getCenterROI(IdealPhase),NUMC,NStep);
    [bestGamma2,loss,exitflag,output] = fminsearch(@myPaperfun,startGamma);
    
    gammaA = 2 * (bestGamma - bestGamma2);
    gammaB = bestGamma - gammaA;
    finalGamma = gammaA / (1 - gammaB);
    endTime = clock;
    
    timecost = etime(endTime,startTime);
end





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
    width = 10;
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

    PECV_first = getPECV(getCenterROI(phaseError), getCenterROI(phaseTure), NUMC, NStep);
    phaseError2 = CompensateC(phaseError,NStep,PECV_first);

    IdealPhase = findPlane(phaseError2);
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
    
    C = A\B;
    
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

function [PhaseImage] = createPhaseImageNoise(Freq,Steps,Width,Height,randomGamma,beta,lambda, m, var)
    PSin = 0; PCos = 0;
    for step = 1:Steps
        n = step - 1;
        
        Im = imnoise(Picture(Freq,Steps,Width,Height,n,randomGamma,beta,lambda),'gaussian',m,var); 
        PSin = PSin + Im * sin(  2 * pi * n / Steps );
        PCos = PCos + Im * cos(  2 * pi * n / Steps );
        clear Im;
    end
    
    PhaseImage = pi + atan2(PSin,  -PCos );
    PhaseImage(:,1) = 0; 
end

function [PhaseImage] = createPhaseImageAmbientLight(Freq,Steps,Width,Height,randomGamma,beta,lambda)
    global alpha_p alpha AmbientLight;

    PSin = 0; PCos = 0;
    for step = 1:Steps
        n = step - 1;
        Im = Picture(Freq,Steps,Width,Height,n,randomGamma,beta,lambda);  % Im [0,1] 
        Im = alpha * alpha_p * Im + AmbientLight;
        PSin = PSin + Im * sin(  2 * pi * n / Steps );
        PCos = PCos + Im * cos(  2 * pi * n / Steps );
        clear Im;
    end
    
    PhaseImage = pi + atan2(PSin,  -PCos );
    PhaseImage(:,1) = 0; 
end

function [PhaseImage] = createPhaseImageDefousing(Freq,Steps,Width,Height,randomGamma,beta,lambda)
    global W;
    PSin = 0; PCos = 0;
    
    for step = 1:Steps
        n = step - 1;
        Im = Picture(Freq,Steps,Width,Height,n,randomGamma,beta,lambda);
        Im = imfilter(Im, W, 'replicate');
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
