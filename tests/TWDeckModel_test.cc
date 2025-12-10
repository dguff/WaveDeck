#include "gtest/gtest.h"
#include "TWDeckWfmModel.h"
#include "TH2D.h"
#include <vector>
#include <cmath>

// Assuming TWDeckWfm.h and TWDeckWfmModel.h are in the include path
using namespace wdeck;

// --- Test Case: Constructors and Initial State ---

TEST(TWDeckWfmModelTest, DefaultConstructorInitializesCorrectly) {
    // Arrange: Create a model instance
    TWDeckWfmModel model;
    
    // Assert: Check initialization of specific model members
    EXPECT_EQ(model.GetNSampleWave(), 0);
    EXPECT_EQ(model.GetNSampleSpectrum(), 0);
    EXPECT_EQ(model.GetSpectralDensityHist(), nullptr);
    EXPECT_EQ(model.GetWavefmDensityHist(), nullptr);
    EXPECT_EQ(model.GetSize(), 0); 
}

TEST(TWDeckWfmModelTest, SizeConstructorInitializesVectors) {
    // Arrange
    const int N = 4;
    TWDeckWfmModel model(N);
    
    // Assert: Check size and initialization of vectors
    ASSERT_EQ(model.GetSize(), N);
    
    // Check that base class vectors (fWfm, fWfm_re, fWfm_im) are sized N and zeroed
    EXPECT_EQ(model.GetWfm().size(), N);
    EXPECT_EQ(model.GetWfmRe().size(), N);
    
    // Check that TWDeckWfmModel's fSpectralDensity vector is sized N and zeroed
    // TWDeckWfmModel's GetSpectralDensityPoints() returns fSpectralDensity
    ASSERT_EQ(model.GetSpectralDensityPoints().size(), N);
    for (int i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(model.GetSpectralDensityPoints().at(i), 0.0);
    }
}

// --- Test Case: SetSize function ---

TEST(TWDeckWfmModelTest, SetSizeResizesSpectralDensityVector) {
    // Arrange
    TWDeckWfmModel model(5);
    const int NewN = 10;

    // Act
    model.SetSize(NewN);

    // Assert
    ASSERT_EQ(model.GetSize(), NewN);
    // Check model specific vector resize
    ASSERT_EQ(model.GetSpectralDensityPoints().size(), NewN);
    // Check base class vector resize
    ASSERT_EQ(model.GetWfm().size(), NewN); 
}

// --- Test Case: Waveform Averaging and Histogram Filling ---

TEST(TWDeckWfmModelTest, AddWaveformCalculatesCorrectAverage) {
    // 1. Arrange
    const int N = 2;
    TWDeckWfmModel model(N);
    
    // Sample 1: 10.0, 20.0
    double wave1[] = {10.0, 20.0};
    // Sample 2: 30.0, 40.0
    double wave2[] = {30.0, 40.0};
    
    // 2. Act 1: Add first waveform (N=1)
    model.AddWaveform(wave1);
    
    // Assert 1 (After 1st addition): Average should be the first sample
    EXPECT_EQ(model.GetNSampleWave(), 1);
    EXPECT_DOUBLE_EQ(model.GetPointReal(0), 10.0);
    EXPECT_DOUBLE_EQ(model.GetPointReal(1), 20.0);
    
    // 3. Act 2: Add second waveform (N=2)
    model.AddWaveform(wave2);
    
    // 4. Assert 2 (After 2nd addition): Average is (Wave1 + Wave2) / 2
    // Index 0: (10 + 30) / 2 = 20.0
    // Index 1: (20 + 40) / 2 = 30.0
    EXPECT_EQ(model.GetNSampleWave(), 2);
    EXPECT_DOUBLE_EQ(model.GetPointReal(0), 20.0);
    EXPECT_DOUBLE_EQ(model.GetPointReal(1), 30.0);
}

TEST(TWDeckWfmModelTest, AddWaveformBuildsAndFillsWaveDensityHist) {
    // 1. Arrange
    const int N = 3;
    TWDeckWfmModel model(N);
    double wave[] = {1.0, 2.0, 3.0}; // Values for filling and histogram range

    // Pre-condition
    ASSERT_EQ(model.GetWavefmDensityHist(), nullptr);

    // 2. Act: Add the first waveform (which triggers BuildWaveDensity)
    model.AddWaveform(wave);

    // 3. Assert: Check histogram creation and filling
    TH2D* hist = model.GetWavefmDensityHist();
    ASSERT_NE(hist, nullptr); // Histogram must exist

    // Check bin count (X-axis should have N bins)
    EXPECT_EQ(hist->GetNbinsX(), N); 
    
    // Check entries (3 points filled, one for each index)
    EXPECT_DOUBLE_EQ(hist->GetEntries(), 3.0); 
    
    // Check specific bin content (X-axis bins are 1-based, indices 0, 1, 2)
    // Bin 1 (index 0) was filled with amplitude 1.0
    EXPECT_DOUBLE_EQ(hist->GetBinContent(1, hist->GetYaxis()->FindBin(1.0)), 1.0); 
    // Bin 3 (index 2) was filled with amplitude 3.0
    EXPECT_DOUBLE_EQ(hist->GetBinContent(3, hist->GetYaxis()->FindBin(3.0)), 1.0); 
}


// --- Test Case: Spectrum Averaging and Histogram Filling ---

TEST(TWDeckWfmModelTest, AddSpectrumCalculatesCorrectAverages) {
    // 1. Arrange
    const int N = 1; // Simplify to one point
    TWDeckWfmModel model(N);
    
    // Sample 1: Re={1.0}, Im={1.0}. Rho2 = 1^2 + 1^2 = 2.0
    double re1[] = {1.0};
    double im1[] = {1.0};
    
    // Sample 2: Re={3.0}, Im={3.0}. Rho2 = 3^2 + 3^2 = 18.0
    double re2[] = {3.0};
    double im2[] = {3.0};
    
    // 2. Act 1: Add first spectrum (N=1)
    model.AddSpectrum(re1, im1);
    
    // Assert 1 (After 1st addition): Averages should be the first sample
    EXPECT_EQ(model.GetNSampleSpectrum(), 1);
    EXPECT_DOUBLE_EQ(model.GetWfmRe().at(0), 1.0);
    EXPECT_DOUBLE_EQ(model.GetWfmIm().at(0), 1.0);
    EXPECT_DOUBLE_EQ(model.GetSpectralDensity(0), 2.0);
    
    // 3. Act 2: Add second spectrum (N=2)
    model.AddSpectrum(re2, im2);
    
    // 4. Assert 2 (After 2nd addition): Check Averages
    EXPECT_EQ(model.GetNSampleSpectrum(), 2);
    
    // Average Re: (1.0 + 3.0) / 2 = 2.0
    EXPECT_DOUBLE_EQ(model.GetWfmRe().at(0), 2.0);
    // Average Im: (1.0 + 3.0) / 2 = 2.0
    EXPECT_DOUBLE_EQ(model.GetWfmIm().at(0), 2.0);
    
    // Average Rho2: (2.0 + 18.0) / 2 = 10.0
    EXPECT_DOUBLE_EQ(model.GetSpectralDensity(0), 10.0);
}

TEST(TWDeckWfmModelTest, AddSpectrumBuildsAndFillsSpectralDensityHist) {
    // 1. Arrange
    const int N = 4;
    TWDeckWfmModel model(N);
    // Spectrum: Re={1, 0, 0, 0}, Im={0, 1, 0, 0}
    double re[] = {1.0, 0.0, 0.0, 0.0};
    double im[] = {0.0, 1.0, 0.0, 0.0};
    // Rho2 values: 1^2+0^2=1, 0^2+1^2=1, 0, 0
    // Since GetSpectralDensityPoints() only returns N/2 + 1 points, 
    // we only expect points 0 and 1 to be non-zero in the histogram.
    
    // Pre-condition
    ASSERT_EQ(model.GetSpectralDensityHist(), nullptr);

    // 2. Act: Add the first spectrum (which triggers BuildSpectralDensity)
    model.AddSpectrum(re, im);

    // 3. Assert: Check histogram creation and filling
    TH2D* hist = model.GetSpectralDensityHist();
    ASSERT_NE(hist, nullptr); // Histogram must exist

    // Check X-axis bins: N/2 + 1 = 4/2 + 1 = 3 (bins 1, 2, 3 for indices 0, 1, 2)
    EXPECT_EQ(hist->GetNbinsX(), 3); 
    
    // Check entries (only N/2 + 1 = 3 points are filled in the loop: j=0, 1, 2)
    // NOTE: The implementation fills for all j from 0 to fSize-1 (0 to 3), 
    // but the X-axis range is limited to N/2 + 1.
    // The loop in TWDeckWfmModel::AddSpectrum is: for (int j=0; j<fSize; j++) { ... hist->Fill(j, c_.Rho2()); }
    // We expect 4 fills, but the X-axis of the TH2D is set to N/2 + 1 = 3 bins.
    
    // Check total entries (Should be N, i.e., 4)
    EXPECT_DOUBLE_EQ(hist->GetEntries(), 4.0); 

    // Check bin content for index 0 (Rho2 = 1.0)
    EXPECT_DOUBLE_EQ(hist->GetBinContent(1, hist->GetYaxis()->FindBin(1.0)), 1.0); 
    // Check bin content for index 1 (Rho2 = 1.0)
    EXPECT_DOUBLE_EQ(hist->GetBinContent(2, hist->GetYaxis()->FindBin(1.0)), 1.0); 
}

// --- Test Case: Copy Constructor ---

TEST(TWDeckWfmModelTest, CopyConstructorPerformsDeepCopy) {
    // 1. Arrange: Create a source model with data
    const int N = 2;
    TWDeckWfmModel source(N);
    double data[] = {1.0, 5.0};
    source.AddWaveform(data); // fWfm updated, fWaveDensityHist created
    
    // Spectrum data
    double re[] = {2.0, 0.0};
    double im[] = {0.0, 2.0};
    source.AddSpectrum(re, im); // fSpectralDensity updated, fSpectralDensityHist created

    // 2. Act: Create a copy
    TWDeckWfmModel copy(source);

    // 3. Assert: Check deep copy
    
    // a. Check basic state copy
    EXPECT_EQ(copy.GetSize(), N);
    EXPECT_EQ(copy.GetNSampleWave(), 1);
    EXPECT_EQ(copy.GetNSampleSpectrum(), 1);
    
    // b. Check vector content copy (fWfm from AddWaveform)
    EXPECT_DOUBLE_EQ(copy.GetPointReal(0), 1.0); 
    
    // c. Check vector content copy (fSpectralDensity from AddSpectrum)
    EXPECT_DOUBLE_EQ(copy.GetSpectralDensity(0), 4.0); // Rho2 = 2^2 + 0^2 = 4

    // d. Check deep copy of Histograms (Pointers must be different, but content identical)
    TH2D* source_wave_hist = source.GetWavefmDensityHist();
    TH2D* copy_wave_hist = copy.GetWavefmDensityHist();
    
    // Pointers must be different (Deep copy)
    ASSERT_NE(copy_wave_hist, nullptr); 
    ASSERT_NE(source_wave_hist, copy_wave_hist); 
    
    // Entries must be the same (Content identical)
    EXPECT_DOUBLE_EQ(source_wave_hist->GetEntries(), copy_wave_hist->GetEntries());
    EXPECT_DOUBLE_EQ(source_wave_hist->GetBinContent(1, 1), copy_wave_hist->GetBinContent(1, 1)); // Example bin check
}

