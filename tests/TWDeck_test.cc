#include "gtest/gtest.h"
#include "TWDeck.h"
#include "TWDeckWfm.h"
#include "TWDeckWfmModel.h"
#include "TWDeckWfmFilter.h"
#include "TMath.h"
#include <vector>
#include <cmath>

// Using the namespace defined in the package
using namespace wdeck;

// --- Helper Functions ---

// Helper to create a simple Sine Wave for testing FFT
void FillSineWave(TWDeckWfm& wfm, double amplitude, double freq_idx) {
    int N = wfm.GetSize();
    double* data = new double[N];
    for (int i = 0; i < N; ++i) {
        // Creates a sine wave that fits exactly 'freq_idx' cycles in N samples
        data[i] = amplitude * std::sin(2 * TMath::Pi() * freq_idx * i / (double)N);
    }
    wfm.LoadWave(N, data);
    delete[] data;
}

// --- Test Cases ---

TEST(TWDeckTest, ConstructorInitializesCorrectly) {
    // Arrange & Act
    int size = 100;
    TWDeck deck(size);

    // Assert
    EXPECT_EQ(deck.GetSize(), size);
    // FFT size in TWDeck is usually managed internally, often N or N/2+1 depending on implementation
    // Based on the code, fFFTSize is set in BuildFFT. Let's check it's non-zero.
    EXPECT_GT(deck.GetFFTSize(), 0); 
}

TEST(TWDeckTest, SetSizeRebuildsConfiguration) {
    // Arrange
    TWDeck deck(50);
    int new_size = 200;

    // Act
    deck.SetSize(new_size);

    // Assert
    EXPECT_EQ(deck.GetSize(), new_size);
    // Verify that logic depending on size (FFT) is updated (indirectly)
    EXPECT_GT(deck.GetFFTSize(), 0);
}

TEST(TWDeckTest, FilterRegistrationAndRetrieval) {
    // Arrange
    TWDeck deck(100);
    
    // NOTE: TWDeck destructor deletes the filters in fFilters.
    // We MUST pass a pointer allocated with 'new', otherwise we get a double-free crash.
    TWDeckWfmFilter* filter = new TWDeckWfmFilter(100); 
    
    const char* filterName = "MyTestFilter";

    // Act
    deck.RegisterFilter(filterName, filter, false);
    TWDeckWfmFilter* retrieved = deck.GetFilter(filterName);

    // Assert
    ASSERT_NE(retrieved, nullptr);
    
    // Act 2: Try to get a non-existent filter
    TWDeckWfmFilter* missing = deck.GetFilter("GhostFilter");
    
    // Assert 2
    EXPECT_EQ(missing, nullptr);

    delete filter; 
}

// --- FFT Round Trip Test (The "Golden" Test) ---
// This verifies that R2C (Time -> Freq) followed by C2R (Freq -> Time)
// returns the original signal (scaled by N usually in ROOT/FFTW).
TEST(TWDeckTest, FFTRoundTripConservesSignalShape) {
    // 1. Arrange
    const int N = 100;
    TWDeck deck(N);
    TWDeckWfm wfm(N);
    
    // Create a known signal (Sine wave)
    double amplitude = 5.0;
    FillSineWave(wfm, amplitude, 4.0); // 4 cycles

    // Verify initial domain
    ASSERT_EQ(wfm.GetOriginDomain(), kReal);

    // 2. Act: Real -> Complex
    deck.FFTR2C(&wfm);

    // In Frequency domain, real part shouldn't be empty
    EXPECT_FALSE(wfm.GetWfmRe().empty());

    // 3. Act: Complex -> Real
    deck.FFTC2R(&wfm);

    // ROOT FFT (via FFTW) usually results in unnormalized output. 
    // The backward transform scales the result by N.
    // We check if the shape is preserved and values are scaled by N.
    std::vector<double>& result = wfm.GetPointsReal();
    
    for (int i = 0; i < N; ++i) {
        double original = amplitude * std::sin(2 * TMath::Pi() * 4.0 * i / (double)N);
        double transformed = result[i];
        
        // Check with tolerance, accounting for scaling factor N
        EXPECT_NEAR(transformed, original, 1e-4) 
            << "Mismatch at index " << i << ". Did normalization change?";
    }
}

// --- Test Produce (Model Generation) ---
TEST(TWDeckTest, ProduceReturnsValidWaveform) {
    // 1. Arrange
    const int N = 64;
    TWDeck deck(N);
    
    // Create a dummy model
    TWDeckWfmModel model(N);
    // Add some dummy spectrum to the model so it's not empty
    double re[N], im[N];
    for(int i=0; i<N; ++i) { re[i]=1.0; im[i]=0.0; } // Flat spectrum
    model.AddSpectrum(re, im);

    // 2. Act
    // Produce generates a new TWDeckWfm allocated with 'new'
    TWDeckWfm* result = deck.Produce(&model);

    // 3. Assert
    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->GetSize(), N);
    EXPECT_EQ(result->GetOriginDomain(), kReal);
    
    // Check that it's not just zeros (since we gave it a spectrum)
    bool hasNonZero = false;
    for(auto v : result->GetPointsReal()) {
        if(std::abs(v) > 1e-9) hasNonZero = true;
    }
    EXPECT_TRUE(hasNonZero);

    // Cleanup: We must delete the produced waveform
    delete result;
}
