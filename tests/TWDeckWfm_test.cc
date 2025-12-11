#include "gtest/gtest.h"
#include "TWDeckWfm.h"
#include <vector>
#include <cmath>

using namespace wdeck;

// --- Test Case: DefaultConstructor ---
TEST(TWDeckWfmTest, DefaultConstructor) {
  TWDeckWfm wfm;

  ASSERT_EQ(wfm.GetSize(), 0);
  ASSERT_EQ(wfm.GetOriginDomain(), kReal);
  ASSERT_TRUE(wfm.GetWfm().empty());
  ASSERT_TRUE(wfm.GetWfmRe().empty());
}

// --- Test Case: SizeConstructor ---
TEST(TWDeckWfmTest, SizeConstructor) {
  // 1. Arrange
  const int N = 10;
  TWDeckWfm wfm(N);

  // 2. Assert: Check size and vector sizes
  ASSERT_EQ(wfm.GetSize(), N);
  ASSERT_EQ(wfm.GetWfm().size(), N);
  ASSERT_EQ(wfm.GetWfmRe().size(), N);

  // Verify that all elements are initialized to 0
  for (int i = 0; i < N; ++i) {
    EXPECT_DOUBLE_EQ(wfm.GetPointReal(i), 0.0);
    EXPECT_DOUBLE_EQ(wfm.GetWfmRe().at(i), 0.0);
  }
}

// --- Test Case: SetSize and ClearWave ---
TEST(TWDeckWfmTest, SetSizeChangesVectors) {
  // 1. Arrange
  TWDeckWfm wfm(5);
  const int NewN = 15;

  // 2. Act: Change size
  wfm.SetSize(NewN);

  // 3. Assert
  ASSERT_EQ(wfm.GetSize(), NewN);
  ASSERT_EQ(wfm.GetWfm().size(), NewN);
  ASSERT_EQ(wfm.GetWfmIm().size(), NewN);

  // Verify that the new elements are zeroed
  for (int i = 0; i < NewN; ++i) {
    EXPECT_DOUBLE_EQ(wfm.GetPointReal(i), 0.0);
  }
}

TEST(TWDeckWfmTest, ClearWaveResetsValues) {
  // 1. Arrange: Load some non-zero data
  const int N = 3;
  TWDeckWfm wfm(N);
  wfm.GetWfm()[0] = 5.0;
  wfm.GetWfmRe()[1] = -2.0;
  wfm.GetWfmIm()[2] = 1.0;

  // 2. Act: Clear the waveform
  wfm.ClearWave();

  // 3. Assert: Check all values are reset to 0
  for (int i = 0; i < N; ++i) {
    EXPECT_DOUBLE_EQ(wfm.GetPointReal(i), 0.0);
    EXPECT_DOUBLE_EQ(wfm.GetWfmRe().at(i), 0.0);
    EXPECT_DOUBLE_EQ(wfm.GetWfmIm().at(i), 0.0);
  }
}

// --- Test Case: LoadWave and LoadSpectrum ---
TEST(TWDeckWfmTest, LoadWaveInitializesRealPart) {
  // 1. Arrange
  const int N = 4;
  double data[] = {1.1, 2.2, 3.3, 4.4};
  TWDeckWfm wfm;

  // 2. Act: Load data array
  wfm.LoadWave(N, data);

  // 3. Assert
  ASSERT_EQ(wfm.GetSize(), N);
  for (int i = 0; i < N; ++i) {
    ASSERT_DOUBLE_EQ(wfm.GetPointReal(i), data[i]);
  }

  // Verify that complex part is zeroed
  for (int i = 0; i < N; ++i) {
    ASSERT_DOUBLE_EQ(wfm.GetWfmRe().at(i), 0.0);
  }
}

TEST(TWDeckWfmTest, LoadSpectrumInitializesComplexPart) {
  // 1. Arrange
  const int N = 2;
  double re_data[] = {10.0, -5.0};
  double im_data[] = {0.5, 9.9};
  TWDeckWfm wfm;

  // 2. Act: Use the method with the size
  wfm.LoadSpectrum(N, re_data, im_data);
  wfm.SetOriginDomain(kComplex); 

  // 3. Assert
  ASSERT_EQ(wfm.GetSize(), N);
  ASSERT_EQ(wfm.GetOriginDomain(), kComplex);

  // Verify the complex coefficients
  TComplex c0 = wfm.GetPointComplex(0);
  TComplex c1 = wfm.GetPointComplex(1);

  ASSERT_DOUBLE_EQ(c0.Re(), 10.0);
  ASSERT_DOUBLE_EQ(c0.Im(), 0.5);
  ASSERT_DOUBLE_EQ(c1.Re(), -5.0);
  ASSERT_DOUBLE_EQ(c1.Im(), 9.9);
}

// --- Test Case: Spectral Density Calculations ---
TEST(TWDeckWfmTest, GetSpectralDensityCalculations) {
  // 1. Arrange
  const int N = 4;
  double re_data[] = {1.0, 3.0, 0.0, -4.0};
  double im_data[] = {2.0, 4.0, 5.0, 0.0};
  // Spectral density (Rho2): Re^2 + Im^2
  // i=0: 1^2 + 2^2 = 5.0
  // i=1: 3^2 + 4^2 = 25.0
  // i=2: 0^2 + 5^2 = 25.0
  // i=3: (-4)^2 + 0^2 = 16.0

  TWDeckWfm wfm;
  wfm.LoadSpectrum(N, re_data, im_data);

  // 2. Assert: Verify individual spectral density values
  EXPECT_DOUBLE_EQ(wfm.GetSpectralDensity(0), 5.0);
  EXPECT_DOUBLE_EQ(wfm.GetSpectralDensity(1), 25.0);
  EXPECT_DOUBLE_EQ(wfm.GetSpectralDensity(2), 25.0);

  // 3. Act & Assert: Test full spectral density vector (size*0.5 + 1 points)
  std::vector<double> density = wfm.GetSpectralDensityPoints();

  // N=4, _size = 4*0.5 + 1 = 3
  ASSERT_EQ(density.size(), 3);
  EXPECT_DOUBLE_EQ(density.at(0), 5.0);
  EXPECT_DOUBLE_EQ(density.at(1), 25.0);
  EXPECT_DOUBLE_EQ(density.at(2), 25.0);
}

