function partV = changeMeshbyEmotion(emotion,partV,partF)

    load('BETA.mat')
    load('sophie.mat')
    load('FacePartsInd.mat')
    load('parsed_data_SL.mat')
    nVertices = 1220;
    disp(size(partV))
    % Unless specified, scale all factors to the same amount 
    scaleFactors = struct('MOUTH', 1, 'EYEL', 1, 'NOSE', 1, 'EYER', 1, 'ALLs', 1);

    if strcmp(emotion, 'smile')
        X = [0.104969700000000;0.104980100000000;0.114777400000000;0.00303261400000000;0.00296590700000000;0.00666228400000000;0.237118100000000;0.238192200000000;-2.18774000000000e-12;-2.18845600000000e-12;0.0560495000000000;0.0560571100000000;0.102377000000000;0.00716774000000000;0;0;0;0;0.685995100000000;0.685998300000000;0.0816218400000000;0.0811506100000000;0.143178100000000;1.53680700000000e-17;0.151617800000000;0.0437056300000000;0.0116943600000000;0.110847700000000;0.110109400000000;0;0;0.0930285600000000;0.0510150900000000;0.792670200000000;0.786838200000000;0.0177960500000000;0.0201759400000000;0.0215139600000000;1.64507500000000e-05;0.0232451400000000;0.0340266800000000;0.00674136400000000;0.287926900000000;0.917157500000000;0.917360100000000;0.212382300000000;0.198032400000000;0.244017900000000;0.245387400000000;0.279677200000000;0.270702900000000;0.000261358700000000];
        X = X*1; 
    elseif strcmp(emotion,'eyes closed')
        X = zeros(52,1); 
        X(9) = 1;
        X(10) = 1;
    elseif strcmp(emotion,'mouth closed')
        X = zeros(52,1); 
        X(27) = 0.05;
    elseif contains(emotion, 'happiness')
        X = sophie.blendshape(:,4472); % displacement  = 5.7642
        if strcmp(emotion, 'happiness_0.5')
            scaleFactors.ALLs = 0.08384; % this value changes manually bc faces are scaled later 
        elseif strcmp(emotion, 'happiness_1')
            scaleFactors.ALLs = 0.16738;
        elseif strcmp(emotion, 'happiness_2')
            scaleFactors.ALLs = 0.33355;
        elseif strcmp(emotion, 'happiness_3')
            scaleFactors.ALLs = 0.49858;
        elseif strcmp(emotion, 'happiness_35')
            scaleFactors.ALLs = 0.58;
        elseif strcmp(emotion, 'happiness_4')
            scaleFactors.ALLs = 0.66262;
        elseif strcmp(emotion, 'happiness_6')
            scaleFactors.ALLs = 0.98832;
        elseif strcmp(emotion, 'happiness_8')
            scaleFactors.ALLs = 1.31205;
        elseif strcmp(emotion, 'happiness_10')
            scaleFactors.ALLs =  1.6336;
        end

    elseif contains(emotion, 'sadness')
        X = sophie.blendshape(:,2173);
        if strcmp(emotion, 'sadness_0.5')
            scaleFactors.ALLs = 0.2745;
        elseif strcmp(emotion, 'sadness_1')
            scaleFactors.ALLs = 0.548;
        elseif strcmp(emotion, 'sadness_2')
            scaleFactors.ALLs = 1.092;
        elseif strcmp(emotion, 'sadness_3')
            scaleFactors.ALLs = 1.632;
        elseif strcmp(emotion, 'sadness_35')
            scaleFactors.ALLs = 1.9;
        elseif strcmp(emotion, 'sadness_4')
            scaleFactors.ALLs = 2.4;
            scaleFactors.MOUTH = 1/1.25;
        % elseif strcmp(emotion, 'sadness_6')
        % 
        %     scaleFactors.ALLs = 2.5;
        %     scaleFactors.EYEL = 4;
        %     scaleFactors.EYER = 4;
        %     scaleFactors.MOUTH = 1/1.2;
        % elseif strcmp(emotion, 'sadness_8')
        %     emotion_scale = 1.39532;
        % elseif strcmp(emotion, 'sadness_10')
        %     emotion_scale =  1.75342;
        end
     

    elseif contains(emotion, 'disgust')
        X = sophie.blendshape(:,3168);
        X(1) = 0.6;
        X(2) = 0.6;
        X(7) = 0.4;
        X(8) = 0.4;
        X(19) = 0.2;
        X(20) = 0.2;
        X(23) = 0.001;
        X(25) = 0.05;
        X(27) = 0.2;
        X(30) = 0.2;
        X(31) = 0.2;
        X(43) = 0.4;
        X(46) = 0.6;
        X(47) = 0.6;
        X(48) = 0.15;
        X(49) = 0.15;
        X(50) = 0.4;
        X(51) = 0.4;
        if strcmp(emotion, 'disgust_0.5')
            scaleFactors.ALLs = 0.227;
        elseif strcmp(emotion, 'disgust_1')
            scaleFactors.ALLs = 0.45;
        elseif strcmp(emotion, 'disgust_2')
            scaleFactors.ALLs = 0.9;
        elseif strcmp(emotion, 'disgust_3')
            scaleFactors.ALLs = 1.49;
            scaleFactors.MOUTH = 1/1.5;
        elseif strcmp(emotion, 'disgust_4')
            scaleFactors.ALLs = 1.9;
            scaleFactors.EYEL = 1.1;
            scaleFactors.EYER = 1.1;
            scaleFactors.NOSE = 1.2;
            scaleFactors.MOUTH = 1/1.45;
        end

    elseif contains(emotion, 'fear')
        X = sophie.blendshape(:,2389);
        X(1) = 0.35;
        X(2) = 0.35;
        X(3) = 0.6;
        X(21) = 0.5;
        X(22) = 0.5;
        X(25) = 0.01;
        X(27) = 0.1;
        X(30) = 0.1;
        X(31) = 0.1;
        X(44) = 0.0;
        X(45) = 0.0;
        X(46) = 0.35;
        X(47) = 0.35;
        if strcmp(emotion, 'fear_0.5')
            scaleFactors.ALLs = 0.30203;
        elseif strcmp(emotion, 'fear_1')
            scaleFactors.ALLs = 0.60285;
        elseif strcmp(emotion, 'fear_2')
            scaleFactors.ALLs = 1.2008;
        elseif strcmp(emotion, 'fear_3')
            scaleFactors.ALLs = 2;
            scaleFactors.EYEL = 1.2;
            scaleFactors.EYER = 1.2;
            scaleFactors.MOUTH = 1/1.5273;
        elseif strcmp(emotion, 'fear_4')
            scaleFactors.ALLs = 3;
            scaleFactors.EYEL = 1.3;
            scaleFactors.EYER = 1.3;
            scaleFactors.MOUTH = 1/1.6;
        end

      

    elseif contains(emotion, 'anger')
        X = sophie.blendshape(:,3181);
        X(1) = 0.6;
        X(2) = 0.6;
        X(3) = 0.0001;
        X(4) = 0.6;
        X(5) = 0.6;
        X(7) = 0;
        X(8) = 0;
        X(19) = 0.;
        X(20) = 0.;
        X(23) = 0;
        X(27) = 0.1;
        X(30) = 0.25;
        X(31) = 0.25;
        X(32) = 0.001;
        X(34) = 0.2;
        X(35) = 0.2;
        X(38) = 0.001;
        X(40) = 0;
        X(42) = 0.7;
        X(43) = 0.5;
        if strcmp(emotion, 'anger_0.5')
            scaleFactors.ALLs = 0.18147;
        elseif strcmp(emotion, 'anger_1')
            scaleFactors.ALLs = 0.35772;
        elseif strcmp(emotion, 'anger_2')
            scaleFactors.ALLs = 0.6911;
        elseif strcmp(emotion, 'anger_3')
            scaleFactors.ALLs = 1.0598;
        elseif strcmp(emotion, 'anger_35')
            scaleFactors.ALLs = 1.3;
        elseif strcmp(emotion, 'anger_4')
            scaleFactors.ALLs = 1.7;
            scaleFactors.MOUTH = 1/1.4;
        end


    elseif contains(emotion, 'surprise')
        X = sophie.blendshape(:,1665);
        X(3) = 0.5;
        X(4) = 0.5;
        X(5) = 0.5;
        X(11) = 0;
        X(12) = 0;
        X(15) = 0.5;
        X(16) = 0.5;
        X(17) = 0.5;
        X(18) = 0.5;
        X(21) = 0.3;
        X(22) = 0.3;
        if strcmp(emotion, 'surprise_0.5')
            scaleFactors.ALLs = 0.05;
            scaleFactors.MOUTH = 1.3072;
        elseif strcmp(emotion, 'surprise_1')
            scaleFactors.ALLs = 0.1;
            scaleFactors.MOUTH = 1.3145;
        elseif strcmp(emotion, 'surprise_2')
            scaleFactors.ALLs = 0.2;
            scaleFactors.MOUTH = 1.3265;
        elseif strcmp(emotion, 'surprise_3')
            scaleFactors.ALLs = 0.3;
            scaleFactors.MOUTH = 1.3389;
        elseif strcmp(emotion, 'surprise_35')
            scaleFactors.ALLs = 0.34;
            scaleFactors.MOUTH = 1.3389;
        elseif strcmp(emotion, 'surprise_4')
            scaleFactors.ALLs = 0.4;
            scaleFactors.MOUTH = 1.3515;
        elseif strcmp(emotion, 'surprise_6')
            scaleFactors.ALLs = 0.72;
            scaleFactors.MOUTH = 1.1038;
        elseif strcmp(emotion, 'surprise_8') % Scale should be 1!
            scaleFactors.ALLs = 1.1625;
            scaleFactors.MOUTH = 1/1.3;
        elseif strcmp(emotion, 'negative_surprise')
            scaleFactors.ALLs = -1;
        end

    elseif strcmp(emotion, 'neutral')
        X = zeros(52,1); 
    end

    X(length(X)+1) = 0;
    delta = reshape(X' * BETA, nVertices, 3);
    
    % Apply scaling factors for each face part
    parts = {'NOSE', 'EYEL', 'EYER', 'MOUTH', 'ALLs'};
    for i = 1:length(parts)
        partName = parts{i};
        partIndices = FacePartsInd.(partName);
        scaleFactor = scaleFactors.(partName);
    
        delta(partIndices, :) = delta(partIndices, :) * scaleFactor;
    end
    
    % Accumulate delta changes
    sub_ind = unique(partF);
    partV = partV + delta(sub_ind, :);
    disp(X);
end