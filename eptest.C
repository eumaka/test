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

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <trackbase_historic/SvtxTrackMap.h>

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

  EpNode += detector;
  if((detector == "EPD") || (detector == "BBC"))
  {
    EpNode += "_SOUTH";
    EventPlaneNodeName.push_back(EpNode);
    EpNode = "EPINFO_" + detector;
    EpNode += "_NORTH";
    EventPlaneNodeName.push_back(EpNode);
  }
  else
  {
    EpNode = "EPINFO_" + detector;
    EventPlaneNodeName.push_back(EpNode);
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

 std::cout << "-------------------------------------------------------------------------------" << std::endl;
 std::cout << "-----------------------Welcome to the Event Plane Module-----------------------" << std::endl;
 std::cout << "-----Please choose from the list: EPD, CEMC, HCALOUT, HCALIN, BBC, TRACKING----" << std::endl;
 std::cout << "The output are the Q vectors and event plane angles up to the default 3rd order" << std::endl;
 std::cout << "-------------------------------------------------------------------------------" << std::endl;
 
 if(detector == "EPD")
 {
 
   std::cout<< "calculating event plane angles for detector " << detector << std::endl;
   TowerNode += detector;
   TowerGeomNode += detector;

   CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
   if (!cent)
   {
      std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
      return;
   }

   int cent_index = cent->get_centile(CentralityInfo::PROP::bimp)/10;
     
   TowerInfoContainerv1 *_towerinfos = findNode::getClass<TowerInfoContainerv1>(topNode, TowerNode);
   if (!_towerinfos)
   {
      std::cout << "Could not locate tower info node " << TowerNode << std::endl;
      exit(1);
   }

   EpdGeom *_epdgeom = findNode::getClass<EpdGeom>(topNode,TowerGeomNode);
   if (!_epdgeom)
   {
      std::cout << "Could not locate geometry node " << TowerGeomNode << std::endl;
      exit(1);
   }


   std::vector<EpHit> epdhitnorth;
   epdhitnorth.clear();

   std::vector<EpHit> epdhitsouth;
   epdhitsouth.clear();
      

   float tile_phi = 0.; float tile_e = 0.; int arm_id = -1;
   int tile_id = -1; int tile_ring = -1; float eMax = 0.;
   std::tuple<unsigned int, unsigned int, unsigned int > tile_index;

   unsigned int ntowers = _towerinfos->size();
   for (unsigned int ch = 0; ch < ntowers;  ch++)
   {

     TowerInfo *_tower = _towerinfos->get_tower_at_channel(ch);
     unsigned int thiskey =_towerinfos->encode_epd(ch);

     tile_phi = _epdgeom->phi(thiskey);
     tile_index = _epdgeom->id_to_side_sector_tile(thiskey);
     arm_id = get<0>(tile_index);
     tile_id = get<2>(tile_index);
     tile_ring = EPDDefs::get_ring(tile_id);

     tile_e = _tower->get_energy();
     eMax = EPDDefs::get_eMax(cent_index, tile_ring);
        
     if(tile_e < 0.2) continue;
     float truncated_tile_e = (tile_e < eMax) ? tile_e : eMax;

     if(arm_id == 0) //south wheel
     {
       EpHit newHit;
       newHit.nMip = truncated_tile_e;
       newHit.phi = tile_phi;
       epdhitsouth.push_back(newHit);
     }
     else if(arm_id == 1) //north wheel
     {
       EpHit newHit;
       newHit.nMip = truncated_tile_e;
       newHit.phi = tile_phi;
       epdhitnorth.push_back(newHit);
     }
          
   }//end loop over sepd tower info
    
 
   EpFinder_det[0]->Results(epdhitsouth, 0, _EpInfo_det[0]);
   EpFinder_det[1]->Results(epdhitnorth, 0, _EpInfo_det[1]);
    
   epdhitsouth.clear();
   epdhitnorth.clear();
   
  }
  
  else if((detector == "CEMC") || (detector == "HCALOUT") || (detector == "HCALIN"))
  {
    
   std::cout<< "calculating event plane angles for detector " << detector << std::endl;
   TowerNode += detector;
   TowerGeomNode += detector;
    
   TowerInfoContainerv1 *_towerinfos = findNode::getClass<TowerInfoContainerv1>(topNode, TowerNode);
   if (!_towerinfos)
   {
     std::cout << "Could not locate tower info node " << TowerNode << std::endl;
     exit(1);
   }

   RawTowerGeomContainer *_towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNode);
   if (!_towergeom)
   {
      std::cout << "Could not locate geometry node " << TowerGeomNode << std::endl;
      exit(1);
   }
 
    
   std::vector<EpHit> calohit;
   calohit.clear();
    
   unsigned int ntowers = _towerinfos->size();
   for (unsigned int ch = 0; ch < ntowers;  ch++)
   {

     TowerInfo *_tower = _towerinfos->get_tower_at_channel(ch);
     unsigned int key =_towerinfos->encode_key(ch);
     RawTowerDefs::CalorimeterId calo_id = RawTowerDefs::convert_name_to_caloid(detector);
     int ieta = _towerinfos->getTowerEtaBin(key);
     int iphi = _towerinfos->getTowerPhiBin(key);
     RawTowerDefs::keytype towerid = RawTowerDefs::encode_towerid(calo_id, ieta, iphi);
     RawTowerGeom *_tgeo = _towergeom->get_tower_geometry(towerid);
     
     if(_tgeo)
     {
       EpHit newHit;
       newHit.nMip = _tower->get_energy();
       newHit.phi = _tgeo->get_phi();
       calohit.push_back(newHit);
     }
  
   }

   EpFinder_det[0]->Results(calohit, 0, _EpInfo_det[0]);
   calohit.clear();
    
  }
  
  else if (detector == "TRACKING")
  {
  
    float px = 0.; float py = 0.; float pz = 0.;
    TVector3 mom;
    
    std::vector<EpHit> trackhit;
    trackhit.clear();
    
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    if(!trackmap)
    {
       std::cout << "Could not locate SvtxTrackMap" << std::endl;
       exit(1);
    }
  
    for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
    {
      SvtxTrack* track = iter->second;
      px = track->get_px();
      py = track->get_py();
      pz = track->get_pz();
      mom.SetX(px); mom.SetY(py); mom.SetZ(pz);
      
      if(mom.Pt() < 1.0)
      {
        EpHit newHit;
        newHit.nMip = 1.;
        newHit.phi = mom.Phi();
        trackhit.push_back(newHit);
      }
    
     }
    
     EpFinder_det[0]->Results(trackhit, 0, _EpInfo_det[0]);
     trackhit.clear();
   }
  
  else if (detector == "BBC")
  {
  
  
  
  
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





