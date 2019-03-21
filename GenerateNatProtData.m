%% Generate sample data for NatProtFig representation

function pts = GenerateNatProtData(FieldSize, RandFieldPoints, PeakSigma, NPeaks, PtsPerPeak)

    randPts = FieldSize*rand(RandFieldPoints, 2);


    peakPts = [];
    peakPost = (FieldSize-2*max(PeakSigma))*rand(NPeaks, 2)+2*max(PeakSigma);
    for k = 1:NPeaks;

        pS = diff(PeakSigma)*rand(1)+min(PeakSigma);
        pN = round(diff(PtsPerPeak)*rand(1)+min(PtsPerPeak));

        peakPts = [peakPts; pS*randn(pN, 2)+repmat(peakPost(k,:), pN, 1)];

    end
    pts = [randPts; peakPts];

end
