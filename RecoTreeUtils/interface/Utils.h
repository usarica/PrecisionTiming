#ifndef UTILS_H
#define UTILS_H

//****************************************************************************************
typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;

//****************************************************************************************
float DeltaPhi(float phi1, float phi2)
{
    float delta_phi = TMath::Abs(phi1 - phi2);
    if (delta_phi > 2*TMath::Pi()) 
	delta_phi -= 2*TMath::Pi();
    if (delta_phi > TMath::Pi() && delta_phi < 2*TMath::Pi()) 
	delta_phi = 2*TMath::Pi() - delta_phi;
    return delta_phi;
}

//****************************************************************************************
float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
    float d_phi = DeltaPhi(phi1, phi2);
    float d_eta = TMath::Abs(eta1 - eta2);

    return TMath::Sqrt(pow(d_eta,2)+pow(d_phi,2));
}

#endif
