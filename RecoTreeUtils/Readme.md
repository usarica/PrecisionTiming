FastTiming
==========

* Tools for event reconstrution exploiting precise timing information, run under CMSSW_6_2_0_SLHC23

* Classes:
  * PFCandidateWithFT: PFCandidate based class, add timing information from SK endcap colorimenter.
  * VertexWithFT: reco::Vertex based class, contains the collection of PFCandidateWithFT related to the vertex.

* Plugins: 
  * RecoFastTiming: 
   - reads the RECO Events tree and produce a ROOT flat trees with the reconstructed time information.
   - the output contains two TTress: one for the particles and one for the vertices.

* Python:
  * The list of samples and file that can be analyzed by RecoFastTiming are listed in RecoFastTiming_cfi.py
  * GetSampleFiles set the right i/o files for the given sample name. Currently available samples:
    * 1) QCD_noPU 
    * 2) QCD_140PU
    * 3) SingleGammaE50_noPU

* In order to run the analyzer execute:
  ```
  cd $CMSSW_BASE/src/FastTiming/RecoTreeUtils/
  scram b -j
  cmsRun test/RecoFastTiming_cfg.py sampleName=QCD_noPU
  ```
