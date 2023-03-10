#include "EpFinderReco.h"

#include "EpFinder.h"
#include "EpInfo.h"    // for EpInfo
#include "EpInfov1.h"  // for EpInfo

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <epd/EPDDefs.h>
#include <epd/EpdGeom.h>
#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>  // for TH2F
#include <TSystem.h>
#include <TVector3.h>  // for TVector3

#include <algorithm>  // for max
#include <cmath>      // for M_PI
#include <cstdlib>    // for NULL, exit, getenv
#include <iostream>
#include <map>  // for _Rb_tree_const_iterator
#include <utility>
#include <vector>  // for vector

using namespace std;

EpFinderReco::EpFinderReco(const std::string &name)
  : SubsysReco(name)
  , detector("EPD")
{
}

EpFinderReco::~EpFinderReco()
{
  for(int i = 0; i < 2; i++)
  {
    delete EpFinder_det[i];
  }
}

int EpFinderReco::Init(PHCompositeNode *topNode)
{
    
  for(int i = 0; i < 2; i++)
  {
    EpFinder_det[i] = new EpFinder(1, 3);
  }
       
    return CreateNodes(topNode);
}

int EpFinderReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *AlgoNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _algonode));
  if (!AlgoNode)
  {
    AlgoNode = new PHCompositeNode(_algonode);
    dstNode->addNode(AlgoNode);
  }

  epnode = "EpInfo_" + detector;
  if( (detector == "EPD") || (detector == "BBC") )
  {
    epnode += "_North";
    EventPlaneNodeName.push_back(epnode);
    epnode = "EpInfo_" + detector;
    epnode += "_South";
    EventPlaneNodeName.push_back(epnode);
  }
  else
  {
    epnode = "EpInfo_" + detector;
    EventPlaneNodeName.push_back(epnode);
  }
    
 for (unsigned int i = 0; i < EventPlaneNodeName.size(); i++)
 {
    EpInfo *EpInfo_det = new EpInfov1();
    PHIODataNode<PHObject> *EpInfo_det_node = new PHIODataNode<PHObject>(EpInfo_det, EventPlaneNodeName[i], "PHObject");
    AlgoNode->addNode(EpInfo_det_node);
 }

 return Fun4AllReturnCodes::EVENT_OK;

}

int EpFinderReco::process_event(PHCompositeNode *topNode)
{

  GetNodes(topNode);
  GetEventPlanes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
    
  for(int i = 0; i < 2; i++)
  {
    EpFinder_det[i]->ResetEvent();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//

void EpFinderReco::GetEventPlanes(PHCompositeNode *topNode)
{

 std::cout << "--------------------------------------------------------------------------------" << std::endl;
 std::cout << "Welcome to Event Plane Module " << std::endl;
 std::cout << "All you need is to choose a detector " << std::endl;
 std::cout << "From the list: EPD, CEMC, HCALOUT, HCALIN, BBC" << std::endl;
 std::cout << "The output are the Q vectors and event plane angles up to the default 3rd order" << std::endl;

 if(detector == "EPD")
 {
 
   std::cout<< "calculating event plane angles for detector " << detector << std::endl;
    
   CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
   if (!cent)
   {
      std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
      return;
   }

   int cent_index = cent->get_centile(CentralityInfo::PROP::bimp)/10;
     
   TowerInfoContainerv1 *_epd_towerinfos_calib = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_EPD");
   if (!_epd_towerinfos_calib)
   {
      std::cout << "Could not locate SEPD CALIB tower info node " << std::endl;
      exit(1);
   }

   EpdGeom *epdtilegeom = findNode::getClass<EpdGeom>(topNode,"TOWERGEOM_EPD");
   if (!epdtilegeom)
   {
      std::cout << "Could not locate SEPD geometry node " << std::endl;
      exit(1);
   }


   std::vector<EpHit> tnehits;
   tnehits.clear();

   std::vector<EpHit> tsehits;
   tsehits.clear();
      

   float tile_phi = 0.; float tile_e = 0.; int arm_id = -1;
   int tile_id = -1; int tile_ring = -1; float eMax = 0.;
   std::tuple<unsigned int, unsigned int, unsigned int > tile_index;

   unsigned int ntowers = _epd_towerinfos_calib->size();
   for (unsigned int ch = 0; ch < ntowers;  ch++)
   {

      TowerInfo *raw_tower = _epd_towerinfos_calib->get_tower_at_channel(ch);
      unsigned int thiskey =_epd_towerinfos_calib->encode_epd(ch);

       tile_phi = epdtilegeom->phi(thiskey);
       tile_index = epdtilegeom->id_to_side_sector_tile(thiskey);
       arm_id = get<0>(tile_index);
       tile_id = get<2>(tile_index);
       tile_ring = EPDDefs::get_ring(tile_id);

       tile_e = raw_tower->get_energy();
       eMax = EPDDefs::get_eMax(cent_index, tile_ring);
        
       if(tile_e < 0.2) continue;
       float truncated_tile_e = (tile_e < eMax) ? tile_e : eMax;

       if(arm_id == 0) //south wheel
     {
            EpHit newHit;
            newHit.nMip = truncated_tile_e;
            newHit.phi = tile_phi;
            tnehits.push_back(newHit);
     }

     else if(arm_id == 1) //north wheel
         {
           EpHit newHit;
           newHit.nMip = truncated_tile_e;
           newHit.phi = tile_phi;
           tsehits.push_back(newHit);
          }
          
      }//end loop over sepd tower info
    
 
     EpFinder_det[0]->Results(tsehits, 0, _EpInfo_det[0]);
     EpFinder_det[1]->Results(tnehits, 0, _EpInfo_det[1]);
    
     tsehits.clear();
     tnehits.clear();
   }
   else if(detector == "CEMC")
   {
     std::cout<<"calculating event plane angles for detector "<< detector<< std::endl;
    
   }
  
  return;
}

int EpFinderReco::GetNodes(PHCompositeNode *topNode)
{
 
   for (unsigned int i = 0; i < EventPlaneNodeName.size(); ++i)
   {
            
    _EpInfo_det[i] = findNode::getClass<EpInfo>(topNode, EventPlaneNodeName[i]);

     if (!_EpInfo_det[i])
     {
      std::cout << PHWHERE << ": Could not find node:"<< EventPlaneNodeName[i] << std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
     }
  }

  return Fun4AllReturnCodes::EVENT_OK;
    
}





