// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANE_EPFINDERRECO_H
#define EVENTPLANE_EPFINDERRECO_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <cmath>

//Forward declarations
class PHCompositeNode;
class EpFinder;
class EpInfo;
class RawTowerGeomContainer;
class PHG4HitContainer;

class EpFinderReco : public SubsysReco
{
 public:
  EpFinderReco(const std::string &name = "EpFinderReco");
  ~EpFinderReco() override;

  int Init(PHCompositeNode *) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  int ResetEvent(PHCompositeNode * /*topNode*/) override;

  void set_algo_node(const std::string &algonode) { _algonode = algonode; }

  void Detector(const std::string &d)
  {
    detector = d;
  }

 private:
 
  void GetEventPlanes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);
  int CreateNodes(PHCompositeNode *);

  std::string _algonode = "EVENT_PLANE";

  std::string detector = "EPD";

  EpFinder *EpFinder_det[2] = {};
  
  EpInfo *_EpInfo_det[2] = {};

  std::vector<std::string> EventPlaneNodeName;
  std::string EpNode = "EPINFO_";
  std::string TowerNode = "TOWERINFO_CALIB_";

};

#endif  //* EVENTPLANE_EPFINDERRECO_H *//
