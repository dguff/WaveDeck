/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck_noise_analysis.cc
 * @created     : Monday Feb 02, 2026 11:11:02 CET
 *
 * \page noise Example of noise analysis and synthesis
 *
 * Source file: `examples/example_wdeck_noise_analysis.cc`
 *
 * This script shows how to use **WaveDeck** to generate synthetic noise waveforms
 * from a reference noise power spectral density (PSD) template through the
 * `TWDeckWfmModel` class.
 *
 * The reference noise PSD template is read from an external text file
 * and used to build a `TWDeckWfmModel` instance. This model is then employed
 * to produce multiple synthetic noise waveforms using the `TWDeck::Produce` method.
 * The generated waveforms are analyzed to compute their average power spectral density,
 * which is then compared to the original template.
 *
 * ![noise analysis example output](example_noise_analysis.png)
 *
 * \include example_wdeck_noise_analysis.cc
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TH1D.h"

#include "TWDeckPATH.h"
#include "TWDeck.h"
#include "TWDeckWfmModel.h"

int example_wdeck_noise_analysis(
    const TString noise_template_path = "../examples/noise_template.txt", 
    const bool debug_waveforms = false)
{
  const bool use_batch_mode = 
    (std::getenv("CI") || std::getenv("GITHUB_ACTIONS")) ? true : false;
  if (use_batch_mode) {
    gROOT->SetBatch(kTRUE);
  }

  // ----------------------------------------------------------
  // setup waveform parameters
  const int wsize = 1024;
  const double sampling_time = 0.016; // [us]
  const double df = 62.5 / 1024; // [MHz]
  const double t0 = 0.0;
  const double t1 = wsize * sampling_time; // [us]

 
  // ----------------------------------------------------------
  // Setup example settings
  const int n_waveforms = 255; 

  // ----------------------------------------------------------
  // Source the reference noise template
  std::ifstream inFile(noise_template_path.Data());
  std::vector<double> noise_template(0.5*wsize, 0.0);
  std::string line_; 
  if (inFile.is_open() == false) {
    std::cerr << "Error opening noise template file!" << std::endl;
    return -1;
  }

  int idx = 0;
  double val_ = 0.0; 
  double freq_ = 0.0; 
  while (std::getline(inFile, line_)) {
    std::stringstream ss; 
    ss << line_;
    ss >> freq_ >> val_;
    noise_template[idx] = val_;
    idx++;
  }
  inFile.close();

  // ----------------------------------------------------------
  // Instance wavedeck
  TWDeck wdeck(wsize);

  // Create noise template model
  TWDeckWfmModel noise_template_model(wsize);
  noise_template_model.SetNameTitle("template_noise_model", "Noise template model");
  noise_template_model.SetOriginDomain(wdeck::kComplex);
  std::vector<double> tmp_re(wsize, 0.0);
  std::vector<double> tmp_im(wsize, 0.0);
  for (int j=0; j<0.5*wsize; j++) {
    tmp_re[j] = sqrt(noise_template[j]);
    tmp_im[j] = 0.0;
    tmp_re[wsize - 1 - j] = sqrt(noise_template[j]);
    tmp_im[wsize - 1 - j] = 0.0;
  }
  noise_template_model.AddSpectrum(&tmp_re[0], &tmp_im[0]);
  
  // ----------------------------------------------------------
  TWDeckWfmModel generated_noise_model(wsize);
  generated_noise_model.SetNameTitle("generated_noise_model", "Generated Noise Sample"); 

  std::vector<TH1D> synthetic_noise_hists = (debug_waveforms) ? 
    std::vector<TH1D>(n_waveforms, 
      TH1D("hSyntheticNoise", "Synthetic Noise Waveform;Time [ns];Amplitude [a.u.]", 
        wsize, t0, t1))
    : std::vector<TH1D>(0);

  for (int iwf = 0; iwf < n_waveforms; ++iwf) {
    // Produce a noise waveform from the template model
    TWDeckWfm* wfm = wdeck.Produce(&noise_template_model);

    if (debug_waveforms == true ) {
      // Store the generated waveform in a histogram for debugging
      auto& h = synthetic_noise_hists.at(iwf);
      h.SetNameTitle(
          Form("hSyntheticNoise_%04d", iwf), 
          Form("Synthetic Noise Waveform %04d;Time [ns];Amplitude [a.u.]", iwf)
          );
      for (int ismp = 0; ismp < wsize; ++ismp) {
        h.SetBinContent(ismp + 1, wfm->GetWfm()[ismp]);
      }
      h.SetLineColor( gStyle->GetColorPalette( iwf * 255 / n_waveforms ) ); 
    }
    
    // Compute the FFT of the generated waveform
    wdeck.FFTR2C(wfm);
    // Add the noise waveform to the model 
    generated_noise_model.AddSpectrum(wfm->GetPointsComplex().data());
    delete wfm;
  }

  gStyle->SetPalette(kViridis);
  gStyle->SetOptStat(0);
  TCanvas* cModel = new TCanvas("cModel", "noise model", 0, 0, 800, 600); 
  cModel->SetTicks(1, 1);
  cModel->SetGrid(1, 1);
  cModel->SetLogy(1);

  auto h2PSD = generated_noise_model.GetSpectralDensityHist(); 
  h2PSD->SetXTitle("Frequency component"); 
  h2PSD->SetYTitle("Power Spectral Density [a.u.]");
  h2PSD->DrawClone("colz"); 

  auto model_pts = noise_template_model.GetSpectralDensityPoints();
  std::vector<double> model_freq(model_pts.size(), 0.0); 
  int ifreq = 0; 
  for (auto& f : model_freq) {f = ifreq; ifreq++;}
  TGraph* gModel = new TGraph(model_pts.size(), model_freq.data(), model_pts.data()); 
  gModel->SetLineWidth(3);
  gModel->DrawClone("l same"); 

  auto avg_pts = generated_noise_model.GetSpectralDensityPoints();
  TGraph* gAvg = new TGraph(avg_pts.size(), model_freq.data(), avg_pts.data());
  gAvg->SetMarkerColor(kGray+1);
  gAvg->SetMarkerStyle(20);
  gAvg->SetMarkerSize(0.5);
  gAvg->DrawClone("p same");

  TLegend* leg = new TLegend(0.11, 0.88, 0.6, 0.75);
  leg->SetMargin(0.07);
  leg->AddEntry(gModel, "Reference Noise Template PSD", "l");
  leg->AddEntry(gAvg, "Generated Noise Avg PSD", "p");
  leg->Draw();
  
  if (debug_waveforms == true) {
    TCanvas* cWfm = new TCanvas("cWfm", "Synthetic Noise Waveforms", 0, 0, 1200, 800);
    for (int iwf = 0; iwf < n_waveforms; ++iwf) {
      auto& h = synthetic_noise_hists.at(iwf);
      TString opt = (iwf == 0) ? "HIST" : "HIST SAME";
      h.DrawClone(opt);
    }
  }

  if (use_batch_mode) {
    cModel->SaveAs("example_4_wdeck_noise_model.png");
  }

  return 0;
}

int main (int argc, char *argv[])
{
  const bool use_batch_mode = 
    (std::getenv("CI") || std::getenv("GITHUB_ACTIONS")) ? true : false;

  TString noise_template_path = Form("%s/examples/noise_generic_model.txt", WDECK_PROJECT_DIR);
  bool debug_waveforms = false;

  if (argc > 1) {
    noise_template_path = TString(argv[1]);
  }
  if (argc > 2) {
    debug_waveforms = (std::atoi(argv[2]) != 0);
  }

  TApplication* tapp = new TApplication("example_wdeck_noise_analysis", &argc, argv);

  example_wdeck_noise_analysis(noise_template_path, debug_waveforms);

  if (!use_batch_mode) tapp->Run();

  return 0;
}

