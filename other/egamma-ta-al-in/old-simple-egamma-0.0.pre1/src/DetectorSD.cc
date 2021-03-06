/* ========================================================== */
// Code was written by Anton Korneev and Mechinsky Vitaly, 
// Institute of Nuclear Problems, Belarus, Minsk, September 2007
//
// And rewritten by Bogdan Maslovskiy,
// Taras Schevchenko National University of Kyiv, 2011.
/* ========================================================== */

#include "DetectorSD.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"

DetectorSD::DetectorSD(G4String name): G4VSensitiveDetector(name)
{
  // получаем указатель на класс RunAction
  // мы будем вызывать его метод RunAction::FillHist
  // для заполнения гистограммы спектра поглощенной энергии
  runAction = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
  
  is_spec_counter = false;
  
  d_hist_min=0;
  d_hist_max=40000;
  d_hist_bins = 1000;
  d_energy_units=1;
  gamma_count=0;
  e_count = 0;
  
  // создаем гистограмму
  // от 0 до 1000, с 1000 каналов
  hist         = new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  photon_hist  = new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  electron_hist= new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  positron_hist= new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  proton_hist  = new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  neutron_hist = new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  alpha_hist   = new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
  other_hist   = new Hist1i(d_hist_min, d_hist_max, d_hist_bins);
}

DetectorSD::~DetectorSD() 
{
  if(hist!=NULL) delete hist;
  if(photon_hist!=NULL) delete photon_hist;
  if(electron_hist!=NULL) delete electron_hist;
  if(positron_hist!=NULL) delete positron_hist;
  if(proton_hist!=NULL) delete proton_hist;
  if(neutron_hist!=NULL) delete neutron_hist;
  if(alpha_hist!=NULL) delete alpha_hist;
  if(other_hist!=NULL) delete other_hist;
}

void DetectorSD::Initialize(G4HCofThisEvent*)
{
  // в начале события сбрасываем энергию поглощенную детектором
  detEnergy = 0;
}

G4bool DetectorSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  if(!is_spec_counter)
    {
      // через детектор летит частица
      // добавляем энергию потерянную частицей
      // к счетчику энергии детектора
      track = step->GetTrack();
      this->particle_name = track->GetDefinition()->GetParticleName();  
      G4double edep = step->GetTotalEnergyDeposit();
      detEnergy += edep;
  
      // G4cout<< "\n---\n"
      // 	<< G4VSensitiveDetector::SensitiveDetectorName
      // 	<< particle_name
      // 	<< " track Deposit energy: "<< edep/keV
      // 	<<" part. kinetic energy:" << track->GetKineticEnergy()/keV
      // 	<<"\n";
    } 
  {
    /**
       otherwise the histrogramm fill functions will be called from
       SteppingAction::UserSteppingAction(const G4Step* aStep),
       if you have previously used
       SteppingAction::SetDetectorSD(std::vector <DetectorSD*> &DSD_vector);
    **/
  }
  return true;
}

void DetectorSD::EndOfEvent(G4HCofThisEvent*)
{
  // сохраняем энергию накопленную за событие в детекторе
  // в гистограмму

  
  //every particles will be written to file "spectrum_total.dat"
  
  //and energy deposit of these particles below will be also 
  //copied to separate file with deposit of particles with given type:
  if (particle_name == "gamma") 
    {fill_hist_photon(detEnergy); return;}
  if (particle_name == "proton") 
    {fill_hist_proton(detEnergy);return;}
  if (particle_name == "neutron") 
    {fill_hist_neutron(detEnergy);return;}
  if (particle_name == "alpha")
    {fill_hist_alpha(detEnergy);return;}
  if (particle_name == "e+") 
    {fill_hist_positron(detEnergy);return;}
  if (particle_name == "e-") 
    {fill_hist_electron(detEnergy);return;}
  
  //other particles will be added to spectrum_total.dat as well
  //without coping to separate file:
  {fill_hist(detEnergy);return;}
    
  
}


/** Set histogram properties.
    \param minimum value of range. 
    Anything lesser than this setpoint will be ignored.
      
    \param maximum value of range. 
    Anything greater than this setpoint will be ignored.
      
    \param Quantity of histogram bins.
*/
void DetectorSD::set_histo(const double min, const double max,
			   const unsigned int n_bins, 
			   const unsigned int E_units)
{
  d_hist_min = min;
  d_hist_max = max;
  if(d_hist_max < d_hist_min)//if somehow on Earth..
    {
      d_hist_max = min;  d_hist_min = max;
    }
  d_hist_bins = n_bins;
  d_energy_units = E_units;
  
  if(hist==NULL) return;
  hist->set(d_hist_min, d_hist_max, d_hist_bins);
  photon_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  electron_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  positron_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  proton_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  neutron_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  alpha_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  other_hist->set(d_hist_min, d_hist_max, d_hist_bins);
  
}

/** private method, it calls fill() method of pointer 
    (*hist_pointer) to object   of class Hist1i.*/
void DetectorSD::fill_hist_pointer(Hist1i *hist_pointer, G4double energy)
{
  if (energy > 0 && (hist_pointer!=NULL))
    {
      switch(d_energy_units)
	{
	case 0:  hist_pointer->fill(energy); break;
	case 1:  hist_pointer->fill(energy/keV); break;
	case 2:  hist_pointer->fill(energy/MeV); break;
	default:
	  hist_pointer->fill(energy/keV); break;
	}
    }
}

/** Get known about particle type from given definition
    and add it's energy to certain histogram
    \param Pointer to particle definition
    \param Energy in keV units.
*/
void DetectorSD::fill_hist(G4ParticleDefinition *pdef, double energy)
{
  if(!pdef || energy < 0) return;
  
  G4String particle_name = pdef->GetParticleName();
  if (particle_name == "gamma") 
    {fill_hist_photon(energy); return;}
  if (particle_name == "proton") 
    {fill_hist_proton(energy);return;}
  if (particle_name == "neutron") 
    {fill_hist_neutron(energy);return;}
  if (particle_name == "alpha")
    {fill_hist_alpha(energy);return;}
  if (particle_name == "e+") 
    {fill_hist_positron(energy);return;}
  if (particle_name == "e-") 
    {fill_hist_electron(energy);return;}
  
  //other particles will be added to spectrum_total.dat as well
  //without coping to separate file:
  {fill_hist(energy);return;}
  
}

/** Get known about particle type from given name
    and add it's energy to certain histogram
    \param Pointer to particle definition
    \param Energy in keV units.
*/
void DetectorSD::fill_hist(const G4String &pname, double energy,unsigned EUNIT)
{
  if(EUNIT>=0 || EUNIT<=2) d_energy_units = EUNIT;
  if(energy < 0) return;
  
  if (pname == "gamma") 
    {fill_hist_photon(energy); return;}
  if (pname == "proton") 
    {fill_hist_proton(energy);return;}
  if (pname == "neutron") 
    {fill_hist_neutron(energy);return;}
  if (pname == "alpha")
    {fill_hist_alpha(energy);return;}
  if (pname == "e+") 
    {fill_hist_positron(energy);return;}
  if (pname == "e-") 
    {fill_hist_electron(energy);return;}
  
  //other particles will be added to spectrum_total.dat as well
  //without coping to separate file:
  {fill_hist(energy);return;}
  
}
  
/** Book particle's energy to the histogram. 
    Calls of methods like 'void fill_hist_electron(double)' will 
    also add values to the same histogram object, that this method 
    is usginh.
    \param particle's energy, eV.
*/
void DetectorSD::fill_hist(G4double energy)
{
  fill_hist_pointer(hist,energy);
}
  
/** Book electron's energy to the histogram.
    \param particle's energy, eV.
*/
void DetectorSD::fill_hist_electron(G4double energy)
{
  fill_hist(energy);
  fill_hist_pointer(electron_hist,energy);
  e_count++;
}

/** Book electron's energy to the histogram.
    \param particle's energy, eV.
*/
void DetectorSD::fill_hist_positron(G4double energy)
{
  fill_hist(energy);
  fill_hist_pointer(positron_hist,energy);
}

/** Book proton's energy to the histogram.
    \param particle's energy, eV.
*/
void DetectorSD::fill_hist_proton(G4double energy)
{
  fill_hist(energy);
  fill_hist_pointer(proton_hist,energy);
}

/** Book alpha particle's energy to the histogram.
    \param particle's energy, eV.
*/
void DetectorSD::fill_hist_alpha(G4double energy)
{
  fill_hist(energy);
  fill_hist_pointer(alpha_hist,energy);
}
  
/** Book neutron's energy to the histogram.
    \param particle's energy, eV.
*/
void DetectorSD::fill_hist_neutron(G4double energy)
{
  fill_hist(energy);
  fill_hist_pointer(neutron_hist,energy);
}
  
/** Book photon's energy to the histogram.
    \param particle's energy
*/
void DetectorSD::fill_hist_photon(G4double energy)
{
  fill_hist(energy);
  fill_hist_pointer(photon_hist,energy);  
  gamma_count++;
}

/** Clear histogramms to start a new data collection.
*/
void DetectorSD::reset_histo()
{
  if(hist!=NULL)         hist       ->clear();
  if(photon_hist!=NULL)  photon_hist->clear();
  if(electron_hist!=NULL)electron_hist->clear();
  if(positron_hist!=NULL)positron_hist->clear();
  if(proton_hist!=NULL)  proton_hist  ->clear();
  if(neutron_hist!=NULL) neutron_hist ->clear();
  if(alpha_hist!=NULL)   alpha_hist->clear();
  if(other_hist!=NULL)   other_hist->clear();
}

/** Write al histogramms to files. 
    It is recommended to call it simultaneously with 
    RunAction::EndOfRunAction(const G4Run* ).
*/
void DetectorSD::save_histo()
{
  // сохраняем гистограмму в файл
  // второй параметр - первая строка файла
  char header_string[256];
  char filename[256];
  
  if(d_energy_units==0)
    sprintf(header_string, "\"energy, eV\", count");
  else
  if(d_energy_units==1)
    sprintf(header_string, "\"energy, %ceV\", count",'k');
  else
  if(d_energy_units==2)
    sprintf(header_string, "\"energy, %ceV\", count",'M');
  
  
  if(hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_total.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      hist->save(filename, header_string);
    }
  if(photon_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_gamma.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      photon_hist->save(filename, header_string);
    }
  if(electron_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_electron.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      electron_hist->save(filename, header_string);
    }
  if(positron_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_positron.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      positron_hist->save(filename, header_string);
    }
  if(proton_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_proton.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      proton_hist->save(filename, header_string);
    }
  if(neutron_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_neutron.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      neutron_hist->save(filename, header_string);
    }

  if(alpha_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_alpha.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      alpha_hist->save(filename, header_string);
    }
      
  if(other_hist!=NULL)
    {
      sprintf(filename, "%s_spectrum_other.hst", G4VSensitiveDetector::SensitiveDetectorName.data());
      other_hist->save(filename, header_string);
    }
  
  G4cout << "\n*********************************************\n";
  G4cout << "\n\n=======================\nEnd of RunAction:\n";
  G4cout << "sensitive reporting: "<< G4VSensitiveDetector::SensitiveDetectorName;
  G4cout << "\nelectrons: "<<e_count;
  G4cout << "\ngamma quants: "<<gamma_count;
  G4cout << "\n=======================\n\n";

}

