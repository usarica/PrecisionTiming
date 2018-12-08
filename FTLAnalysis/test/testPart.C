{
edm::Wrapper<std::vector<reco::GenParticle>>* wrapperptr=nullptr;
Events->SetBranchStatus("recoGenParticles*", 1); Events->SetBranchAddress("recoGenParticles_genParticles__HLT.", &wrapperptr);
Events->GetEntry(0);


for (int ev=0; ev<Events->GetEntries(); ev++){
  wrapperptr=nullptr;
  Events->ResetBranchAddresses(); Events->SetBranchAddress("recoGenParticles_genParticles__HLT.", &wrapperptr);
  Events->GetEntry(ev);
  std::cout << "Processing entry " << ev << " / " << Events->GetEntries() << std::endl;
  if (!wrapperptr || !wrapperptr->isPresent()){ std::cout << "\t- Gen info not found!" << std::endl; continue; }
  else{ std::cout << "There are " << wrapperptr->product()->size() << " gen. parts." << std::endl; }

  for (unsigned int ip=0; ip<wrapperptr->product()->size(); ip++){
    if (
      abs(wrapperptr->product()->at(ip).pdgId())==13
      &&
      !wrapperptr->product()->at(ip).isPromptFinalState()
      &&
      wrapperptr->product()->at(ip).numberOfMothers()>0
      &&
      wrapperptr->product()->at(ip).mother(0)
      &&
      abs(wrapperptr->product()->at(ip).mother(0)->pdgId())==24
      &&
      wrapperptr->product()->at(ip).mother(0)->numberOfMothers()>0
      &&
      wrapperptr->product()->at(ip).mother(0)->mother(0)
      ) std::cout << wrapperptr->product()->at(ip).mother(0)->mother(0)->pdgId() << std::endl;
  }
}
}