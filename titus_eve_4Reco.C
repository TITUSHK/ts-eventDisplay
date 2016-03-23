// Based on alice_esd.C from ROOT tutorials
// Modified for titus by Alex Finch

#include <iostream>

#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TRootBrowser.h"
#include "TStyle.h"

#include "TGFileDialog.h"
#include "TGButton.h"
#include "TGTab.h"

#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoTube.h"
#include "TGeoShape.h"
#include "TGeoSphere.h"

#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TEveElement.h"
#include "TEveEventManager.h"
#include "TEveScene.h"
#include "TEveSceneInfo.h"
#include "TEveViewer.h"
#include "TEveTrack.h"
#include "TEveTrans.h"
#include "TEveBoxSet.h"
#include "TEveBox.h"
#include "TEveGeoShape.h"
#include "TEveGeoNode.h"
#include "TEveRGBAPalette.h"
#include "TTimer.h"

#include "TGLViewer.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include <map>


using namespace std;


const float tankR=5500;
const float tankZ=11000;

bool haveReconFile;
bool drawPrimaryOnly;
Double_t maxX,maxY,maxZ,minZ;

TTree * TitusEvents;
TTree * TitusPMTs;
TTree * Low_E;
TTree * High_E_Electron;
TTree * High_E_Muon;
TTree * Final_Reconstruction;

void       make_gui();
void       load_event();
void       loadReconstructedEvent();
TGCheckButton *fColourIsTimeBox;
TGeoVolume *WorldVolume,*SimpleVolume;
TEveElementList *FlatGeometry;
TEveScene*  UnrolledScene;
TEveScene*  flatGeometryScene;
TEveViewer* UnrolledView;
TGNumberEntry* 		NumberEntry;
TGNumberEntry* 		NEUTMode;
TGLabel* NEUTModeLabel;
TEveRGBAPalette *pal3D,*palUnrolled;
TEveViewer* New2dView(TString name,TGLViewer::ECameraType type, TEveScene* scene);
bool loadPMT(int hit_PMTid);
void     CreateCube(float vert[24],TLorentzVector centreVec,float side);
void               NextCubeVertex(float vert[6],float step,int changeMe,float sign);

class PaletteHandler;
std::vector<PaletteHandler*> allPaletteHandlers;
Int_t event_id       = 0; // Current event id.

TEveTrackList *gTrackList = 0;

// Declaration of leaf types
Int_t           evt;
#define MAXHITS 1000000

Int_t           nhits;
Double_t        hit_time[MAXHITS];   //[nhits]
Double_t        hit_x[MAXHITS];   //[nhits]
Double_t        hit_y[MAXHITS];   //[nhits]
Double_t        hit_z[MAXHITS];   //[nhits]
Int_t           hit_PMTid[MAXHITS];   //[nhits]

#define MAXPART 1000000
Int_t           npart;
Double_t        part_xStart[MAXPART];   //[npart]
Double_t        part_yStart[MAXPART];   //[npart]
Double_t        part_zStart[MAXPART];   //[npart]
Double_t        part_tStart[MAXPART];   //[npart]
Double_t        part_xEnd[MAXPART];   //[npart]
Double_t        part_yEnd[MAXPART];   //[npart]
Double_t        part_zEnd[MAXPART];   //[npart]
Double_t        part_tEnd[MAXPART];   //[npart]
Double_t        part_pxStart[MAXPART];   //[npart]
Double_t        part_pyStart[MAXPART];   //[npart]
Double_t        part_pzStart[MAXPART];   //[npart]
Double_t        part_pxEnd[MAXPART];   //[npart]
Double_t        part_pyEnd[MAXPART];   //[npart]
Double_t        part_pzEnd[MAXPART];   //[npart]
Double_t        part_KEstart[MAXPART];   //[npart]
Double_t        part_KEend[MAXPART];   //[npart]
Int_t           part_processStart[MAXPART];   //[npart]
Int_t           part_processEnd[MAXPART];   //[npart]
Int_t           part_parentid[MAXPART];   //[npart]
Int_t           part_trackid[MAXPART];   //[npart]
Int_t           part_pid[MAXPART];   //[npart]

Int_t           ncapturecount;
Int_t           neutroncount;
#define  MAXCAPTURES 500
Double_t        capt_x[MAXCAPTURES];   //[ncapturecount]
Double_t        capt_y[MAXCAPTURES];   //[ncapturecount]
Double_t        capt_z[MAXCAPTURES];   //[ncapturecount]
Double_t        capt_t0[MAXCAPTURES];   //[ncapturecount]
Double_t        capt_E[MAXCAPTURES];   //[ncapturecount]
Int_t           capt_num[MAXCAPTURES];   //[ncapturecount]
Int_t           capt_pid[MAXCAPTURES];   //[ncapturecount]
Int_t           capt_nucleus[MAXCAPTURES];   //[ncapturecount]
Int_t           capt_nphot[MAXCAPTURES];   //[ncapturecount]
Int_t           capt_ngamma[MAXCAPTURES];   //[ncapturecount]

Int_t           mode;
Double_t        neutrino_E;
Int_t           neutrino_id;
Double_t        neutrino_px;
Double_t        neutrino_py;
Double_t        neutrino_pz;
Int_t           ntrks;
Int_t           nneutrons;
Double_t        vtxx;
Double_t        vtxy;
Double_t        vtxz;

#define MAXTRKS  100
Int_t           mpid[MAXTRKS];   //[ntrks]
Double_t        px[MAXTRKS];   //[ntrks]
Double_t        py[MAXTRKS];   //[ntrks]
Double_t        pz[MAXTRKS];   //[ntrks]
Double_t        KE[MAXTRKS];   //[ntrks]

// List of branches
TBranch        *b_evt;   //!
TBranch        *b_nhits;   //!
TBranch        *b_hit_time;   //!
TBranch        *b_hit_x;   //!
TBranch        *b_hit_y;   //!
TBranch        *b_hit_z;   //!
TBranch        *b_hit_PMTid;   //!
TBranch        *b_npart;   //!
TBranch        *b_part_xStart;   //!
TBranch        *b_part_yStart;   //!
TBranch        *b_part_zStart;   //!
TBranch        *b_part_tStart;   //!
TBranch        *b_part_xEnd;   //!
TBranch        *b_part_yEnd;   //!
TBranch        *b_part_zEnd;   //!
TBranch        *b_part_tEnd;   //!
TBranch        *b_part_pxStart;   //!
TBranch        *b_part_pyStart;   //!
TBranch        *b_part_pzStart;   //!
TBranch        *b_part_pxEnd;   //!
TBranch        *b_part_pyEnd;   //!
TBranch        *b_part_pzEnd;   //!
TBranch        *b_part_KEstart;   //!
TBranch        *b_part_KEend;   //!
TBranch        *b_part_processStart;   //!
TBranch        *b_part_processEnd;   //!
TBranch        *b_part_parentid;   //!
TBranch        *b_part_trackid;   //!
TBranch        *b_part_pid;   //!
TBranch        *b_ncapturecount;   //!
TBranch        *b_neutroncount;   //!
TBranch        *b_capt_x;   //!
TBranch        *b_capt_y;   //!
TBranch        *b_capt_z;   //!
TBranch        *b_capt_t0;   //!
TBranch        *b_capt_E;   //!
TBranch        *b_capt_num;   //!
TBranch        *b_capt_pid;   //!
TBranch        *b_capt_nucleus;   //!
TBranch        *b_capt_nphot;   //!
TBranch        *b_capt_ngamma;   //!
TBranch        *b_mode;   //!
TBranch        *b_neutrino_E;   //!
TBranch        *b_neutrino_id;   //!
TBranch        *b_neutrino_px;   //!
TBranch        *b_neutrino_py;   //!
TBranch        *b_neutrino_pz;   //!
TBranch        *b_ntrks;   //!
TBranch        *b_nneutrons;   //!
TBranch        *b_vtxx;   //!
TBranch        *b_vtxy;   //!
TBranch        *b_vtxz;   //!
TBranch        *b_mpid;   //!
TBranch        *b_px;   //!
TBranch        *b_py;   //!
TBranch        *b_pz;   //!
TBranch        *b_KE;   //!


/*
Now PMTsTre
*/
// Declaration of leaf types
Int_t           pmt_id;
Double_t        pmt_pos_x;
Double_t        pmt_pos_y;
Double_t        pmt_pos_z;
Double_t        pmt_size;
Double_t        pmt_qe;
Double_t        pmt_time_res;
Int_t           is_lappd;


// List of branches
TBranch        *b_pmt_id;   //!
TBranch        *b_pmt_pos_x;   //!
TBranch        *b_pmt_pos_y;   //!
TBranch        *b_pmt_pos_z;   //!
TBranch        *b_pmt_size;   //!
TBranch        *b_pmt_qe;   //!
TBranch        *b_pmt_time_res;   //!
TBranch        *b_is_lappd;   //!

/*
now the reconstructed information, first Low_E
*/
// Declaration of leaf types
Int_t           Low_E_evt;
Int_t           Low_E_nClusters;
Int_t           Low_E_nSubevents;
Int_t           Low_E_cluster[20];   //[nClusters]
Double_t        recoVtxXLowE[20];   //[nClusters]
Double_t        recoVtxYLowE[20];   //[nClusters]
Double_t        recoVtxZLowE[20];   //[nClusters]
Double_t        recoTimeLowE[20];   //[nClusters]
Double_t        recoDirXLowE[20];   //[nClusters]
Double_t        recoDirYLowE[20];   //[nClusters]
Double_t        recoDirZLowE[20];   //[nClusters]
Double_t        recoChkvAngleLowE[20];   //[nClusters]
Double_t        recoEnergyLowE[20];   //[nClusters]

// List of branches
TBranch        *b_Low_E_evt;   //!
TBranch        *b_Low_E_nClusters;   //!
TBranch        *b_Low_E_nSubevents;   //!
TBranch        *b_Low_E_cluster;   //!
TBranch        *b_recoVtxXLowE;   //!
TBranch        *b_recoVtxYLowE;   //!
TBranch        *b_recoVtxZLowE;   //!
TBranch        *b_recoTimeLowE;   //!
TBranch        *b_recoDirXLowE;   //!
TBranch        *b_recoDirYLowE;   //!
TBranch        *b_recoDirZLowE;   //!
TBranch        *b_recoChkvAngleLowE;   //!
TBranch        *b_recoEnergyLowE;   //!

/*
High_E_Electron
*/
// Declaration of leaf types
Int_t           High_E_Electron_evt;
Int_t           High_E_Electron_nClusters;
Int_t           High_E_Electron_nSubevents;
Int_t           High_E_Electron_cluster[20];   //[nSubevents]
Int_t           High_E_Electron_ring[20];   //[nSubevents]
Double_t        recoVtxXHighEElectron[20];   //[nSubevents]
Double_t        recoVtxYHighEElectron[20];   //[nSubevents]
Double_t        recoVtxZHighEElectron[20];   //[nSubevents]
Double_t        recoTimeHighEElectron[20];   //[nSubevents]
Double_t        recoDirXHighEElectron[20];   //[nSubevents]
Double_t        recoDirYHighEElectron[20];   //[nSubevents]
Double_t        recoDirZHighEElectron[20];   //[nSubevents]
Double_t        recoChkvAngleHighEElectron[20];   //[nSubevents]
Double_t        recoEnergyHighEElectron[20];   //[nSubevents]
Double_t        recoLnLHighEElectron[20];   //[nSubevents]

// List of branches
TBranch        *b_High_E_Electron_evt;   //!
TBranch        *b_High_E_Electron_nClusters;   //!
TBranch        *b_High_E_Electron_nSubevents;   //!
TBranch        *b_High_E_Electron_cluster;   //!
TBranch        *b_High_E_Electron_ring;   //!
TBranch        *b_recoVtxXHighEElectron;   //!
TBranch        *b_recoVtxYHighEElectron;   //!
TBranch        *b_recoVtxZHighEElectron;   //!
TBranch        *b_recoTimeHighEElectron;   //!
TBranch        *b_recoDirXHighEElectron;   //!
TBranch        *b_recoDirYHighEElectron;   //!
TBranch        *b_recoDirZHighEElectron;   //!
TBranch        *b_recoChkvAngleHighEElectron;   //!
TBranch        *b_recoEnergyHighEElectron;   //!
TBranch        *b_recoLnLHighEElectron;   //!

/*
High E muon
*/
// Declaration of leaf types
Int_t           High_E_Muon_evt;
Int_t           High_E_Muon_nClusters;
Int_t           High_E_Muon_nSubevents;
Int_t           High_E_Muon_cluster[20];   //[nSubevents]
Int_t           High_E_Muon_ring[20];   //[nSubevents]
Double_t        recoVtxXHighEMuon[20];   //[nSubevents]
Double_t        recoVtxYHighEMuon[20];   //[nSubevents]
Double_t        recoVtxZHighEMuon[20];   //[nSubevents]
Double_t        recoTimeHighEMuon[20];   //[nSubevents]
Double_t        recoDirXHighEMuon[20];   //[nSubevents]
Double_t        recoDirYHighEMuon[20];   //[nSubevents]
Double_t        recoDirZHighEMuon[20];   //[nSubevents]
Double_t        recoChkvAngleHighEMuon[20];   //[nSubevents]
Double_t        recoEnergyHighEMuon[20];   //[nSubevents]
Double_t        recoLnLHighEMuon[20];   //[nSubevents]

// List of branches
TBranch        *b_High_E_Muon_evt;   //!
TBranch        *b_High_E_Muon_nClusters;   //!
TBranch        *b_High_E_Muon_nSubevents;   //!
TBranch        *b_High_E_Muon_cluster;   //!
TBranch        *b_High_E_Muon_ring;   //!
TBranch        *b_recoVtxXHighEMuon;   //!
TBranch        *b_recoVtxYHighEMuon;   //!
TBranch        *b_recoVtxZHighEMuon;   //!
TBranch        *b_recoTimeHighEMuon;   //!
TBranch        *b_recoDirXHighEMuon;   //!
TBranch        *b_recoDirYHighEMuon;   //!
TBranch        *b_recoDirZHighEMuon;   //!
TBranch        *b_recoChkvAngleHighEMuon;   //!
TBranch        *b_recoEnergyHighEMuon;   //!
TBranch        *b_recoLnLHighEMuon;   //!
/*
Final Reconstruction
*/

// Declaration of leaf types
Int_t           Final_evt;
Int_t           Final_nClusters;
Int_t           Final_nSubevents;
Int_t           recoNRings[20];   //[nClusters]
Int_t           Final_cluster[20];   //[nSubevents]
Int_t           Final_ring[20];   //[nSubevents]
Double_t        recoVtxX[20];   //[nSubevents]
Double_t        recoVtxY[20];   //[nSubevents]
Double_t        recoVtxZ[20];   //[nSubevents]
Double_t        recoTime[20];   //[nSubevents]
Double_t        recoDirX[20];   //[nSubevents]
Double_t        recoDirY[20];   //[nSubevents]
Double_t        recoDirZ[20];   //[nSubevents]
Double_t        recoChkvAngle[20];   //[nSubevents]
Double_t        recoEnergy[20];   //[nSubevents]
Int_t           recoPID[20];   //[nSubevents]

// List of branches
TBranch        *b_Final_evt;   //!
TBranch        *b_Final_nClusters;   //!
TBranch        *b_Final_nSubevents;   //!
TBranch        *b_recoNRings;   //!
TBranch        *b_Final_cluster;   //!
TBranch        *b_Final_ring;   //!
TBranch        *b_recoVtxX;   //!
TBranch        *b_recoVtxY;   //!
TBranch        *b_recoVtxZ;   //!
TBranch        *b_recoTime;   //!
TBranch        *b_recoDirX;   //!
TBranch        *b_recoDirY;   //!
TBranch        *b_recoDirZ;   //!
TBranch        *b_recoChkvAngle;   //!
TBranch        *b_recoEnergy;   //!
TBranch        *b_recoPID;   //!

struct hitStore
{
	int count;
	float time;
};
struct hitStoreCompare{
	bool operator()(pair<int, hitStore> a,pair<int, hitStore> b)
	{
		return a.second.time<b.second.time;
	}
};	
void       UnrollView(double* pmtX ,double* pmtY,double* pmtZ,int location,float maxY,float maxZ);
int PMTLocation()
{
	//if(R<5500)
	//{
	//	cout<<" tube is at "<<pmt_pos_x<<" "<<pmt_pos_y<<" "<<pmt_pos_z<<endl;
	//	cout<<" so R = "<<sqrt(pmt_pos_x*pmt_pos_x+pmt_pos_y*pmt_pos_y)<<endl;
	//}
	int  location=999;
	if(TMath::Abs(pmt_pos_z-maxZ)<0.1){
		location=0;
	}
	else if(TMath::Abs(pmt_pos_z+maxZ)<0.1){
		location=2;
	}
	else{
		location=1;
	}
	return location;
}
class PaletteHandler
{
	TTimer*              fTimer;
	int  minValue,maxValue,  lowLimit,highLimit;
	TString me;
	TEveRGBAPalette* fPalette;
	unsigned int myNumber;
	static bool busy;
public:
	PaletteHandler(TString m,TEveRGBAPalette* f) : me(m), fPalette(f)  
	{
		fTimer = new TTimer();
		fTimer->Connect("Timeout()", "PaletteHandler", this,"TimedOut()");
		myNumber=allPaletteHandlers.size();
		allPaletteHandlers.push_back(this);		
		busy=kFALSE;
	}
	void setMinMax(int minValue2Set,int maxValue2Set)
	{
		fPalette->SetMinMax(minValue2Set,maxValue2Set);
		minValue=minValue2Set;
		maxValue=maxValue2Set;
	}
	TString Called(){return me;}
	void TimedOut()
	{
		if(busy)return;
		busy=kTRUE;
		for(unsigned int i=0;i<allPaletteHandlers.size();i++)
		{
			if(i!=myNumber){
				allPaletteHandlers[i]->setMinMax(minValue,maxValue);
			}
		}
		gEve->GetScenes()->RepaintAllScenes(kFALSE);
		gEve->DoRedraw3D();
		busy=kFALSE;
	}
	void PaletteChanged()
	{
		minValue=fPalette->GetMinVal();
    maxValue=fPalette->GetMaxVal();
    fTimer->Start(500,kTRUE);
	}
};
bool PaletteHandler::busy;

/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/
void createGeometry(bool flatTube)
{
	
	/*
	Scan the list of PMTs to find the range in x,y,z
	*/
	maxX=0;
	maxY=0;
	maxZ=0;
	minZ=0; 
	for(int pmt=0;pmt<TitusPMTs->GetEntries();pmt++)
	{
		//		cout<<pmt<<" "<<pmt_id<<" "<<pmt_pos_x<<" "<<pmt_pos_y<<" "<<pmt_pos_z<<endl;
		TitusPMTs->GetEntry(pmt);
		maxX=TMath::Max(maxX,pmt_pos_x);
		maxY=TMath::Max(maxY,pmt_pos_y);
		maxZ=TMath::Max(maxZ,pmt_pos_z);
	}
	gSystem->Load("libGeom");
	new TGeoManager("Titus", "Titus Detector");
	// MATERIALS, MIXTURES AND TRACKING MEDIA
	// Material: world
	Double_t a       = 0.000000;
	Double_t z       = 0.000000;
	Double_t density = 0.000000;
	Double_t radl    = 1000000000000000019884624838656.000000;
	Double_t absl    = 1000000000000000019884624838656.000000;
	TGeoMaterial* pMat1 = new TGeoMaterial("world", a,z,density,radl,absl);
	pMat1->SetIndex(0);
	pMat1->SetTransparency(60);
	// Medium: medium0
	Int_t numed   = 0;  // medium number
	Double_t par[20];
	par[0]  = 0.000000; // isvol
	par[1]  = 0.000000; // ifield
	par[2]  = 0.000000; // fieldm
	par[3]  = 0.000000; // tmaxfd
	par[4]  = 0.000000; // stemax
	par[5]  = 0.000000; // deemax
	par[6]  = 0.000000; // epsil
	par[7]  = 0.000000; // stmin
	TGeoMedium* pMed1 = new TGeoMedium("medium0", numed,pMat1, par);
	
	
	float dx = 2*maxX;
	float dy = 2*maxY;
	float dz = 2*maxZ;
	
	
	TGeoShape *pworldbox_1 = new TGeoBBox("worldbox", dx,dy,dz);
	// Volume: volume0
	WorldVolume = new TGeoVolume("FullGeometry",pworldbox_1, pMed1);
	WorldVolume->SetVisLeaves(kTRUE);
	// Volume: volume0
	SimpleVolume = new TGeoVolume("SimpleGeometry",pworldbox_1, pMed1);
	SimpleVolume->SetVisLeaves(kTRUE);
	
	
	
	FlatGeometry = new TEveElementList("Flat Geometry");
	
	TEveGeoShape *Shape;
	
	gGeoManager->SetMaxVisNodes(TitusPMTs->GetEntries());
	if(flatTube)	
		gGeoManager->SetTopVolume(SimpleVolume);
	else	
		gGeoManager->SetTopVolume(WorldVolume);
	
	gGeoManager->GetTopNode()->ls();
	float PMTRadius = pmt_size/2.0;
	cout<<" Phototube radius is "<<PMTRadius*0.0393700787<<" inches "<<endl;
	Double_t rmin   = 0.900000*PMTRadius;
	Double_t rmax   = PMTRadius;
	Double_t theta1 = 0.000000;
	Double_t theta2 = 90.000000;
	Double_t phi1   = 0.000000;
	Double_t phi2   = 360.000000;
	/*
	Full geometry uses hemisphere for tube,  flatTube means just a disk
	*/
	Double_t zSize=0.1*PMTRadius;
	TGeoShape *PhotoTubeShape   = new TGeoTube("phototube",0.0,rmax,zSize);  
	TGeoShape *PhotoSphereShape = new TGeoSphere("photosphere",rmin,rmax,theta1, theta2,phi1,phi2);
	// Volume: phototube
	TGeoVolume* PhotoTubeVolume = new TGeoVolume("phototube",PhotoTubeShape, pMed1);
	PhotoTubeVolume->SetVisLeaves(kTRUE);
	PhotoTubeVolume->SetLineColor(kYellow-5);
	
	TGeoVolume* PhotoSphereVolume = new TGeoVolume("photosphere",PhotoSphereShape, pMed1);
	PhotoSphereVolume->SetVisLeaves(kTRUE);
	PhotoSphereVolume->SetLineColor(kYellow-5);
	TGeoRotation rota("rot",10,20,30);
	TGeoTranslation trans;
	//TGeoCombiTrans *c1 = new TGeoCombiTrans(trans,rota);
	/*
	Read in location of ALL phototubes
	*/
	Double_t theta, phi;
	double pmtX2=0;
	double pmtY2=0;
	double pmtZ2=0;
	//float oldx=0;
	//float oldy=0;
	//float oldz=0;
	for(int tubeId=0;tubeId<TitusPMTs->GetEntries();tubeId++)
	{
		TitusPMTs->GetEntry(tubeId);
		
		//cout<<tubeId<<" "<<pmt_pos_x<<" "<<pmt_pos_y<<" "<<pmt_pos_z<<" "<<pmt_size;
		//if(tubeId!=0)
		//{
	 	//	float area=3.14159*(pmt_size/2)*(pmt_size/2);
		//	float step=sqrt((oldx-pmt_pos_x)*(oldx-pmt_pos_x)+(oldy-pmt_pos_y)*(oldy-pmt_pos_y)+(oldz-pmt_pos_z)*(oldz-pmt_pos_z));
		//	std::cout<<" "<<step<<" ";
		//  std::cout<<area<<" "<<100*area/(step*step);
		//}
		//  oldx=pmt_pos_x;
		//  oldy=pmt_pos_y;
		//  oldz=pmt_pos_z;
		//std::cout<<std::endl;
		double pmtX = pmt_pos_x;
		double pmtY = pmt_pos_y;
		double pmtZ = pmt_pos_z;
		pmtX2=pmtX;
		pmtY2=pmtY;
		pmtZ2=pmtZ;
		int location=PMTLocation();
		UnrollView(&pmtX2,&pmtY2,&pmtZ2,location,maxY,maxZ);
		
		theta=0.0;
		phi=0.0;
		float rad2deg=57.3;
		if(location==1)
		{
			float lengthXY=sqrt(pmtX*pmtX+pmtY*pmtY);
			float xd=pmtX/lengthXY;
			float yd=pmtY/lengthXY;
			TVector3 D(xd,yd,0.0);
			theta=D.Theta()*rad2deg;
			phi=D.Phi()*rad2deg;
			//	cout<<" barrel tube so theta = "<<theta<<" and phi = "<<phi<<endl;
		}
		else
		{
			theta=0.0;
			phi=0.0;
			if(location==2)theta=180.0; 
		}
		//cout<<" tube> "<<tubeId<<endl;
		if(location==0 || location ==1 || location ==2){
			/* create a fake geometry object out of tevegeoshapes for the rolled out view*/
			//cout<<" create translation using "<<pmtX2<<" "<<pmtY2<<" "<<pmtZ2<<endl;
			TGeoTranslation PhototubeUnrolledPositionMatrix("ExlodedShift",pmtX2,pmtY2,pmtZ2);
			Shape = new TEveGeoShape(Form("Phototube %i",tubeId));
			Shape->SetShape(PhotoTubeShape);
			Shape->SetTransMatrix(PhototubeUnrolledPositionMatrix);
			Shape->SetMainColor(kYellow-5);
			Shape->SetMainTransparency(70);
			FlatGeometry->AddElement(Shape);
			/* now the 'normal' root geometry objects */
			TGeoRotation TubeRotation("rotation",phi-90.0,theta,0.0);//D.Phi(),D.Theta,0.0);
			TGeoTranslation PhototubePositionMatrix("shift",pmtX,pmtY,pmtZ);
			TGeoCombiTrans *ShiftAndTwist = new TGeoCombiTrans(PhototubePositionMatrix,TubeRotation);
			SimpleVolume->AddNode(PhotoTubeVolume, tubeId, ShiftAndTwist);
			WorldVolume-> AddNode(PhotoSphereVolume, tubeId, ShiftAndTwist);
		}
		
	}
	// CLOSE GEOMETRY
	gGeoManager->CloseGeometry();
	
}
// Convert x,y,z of phototubes to position in Unrolled view
void       UnrollView(double* pmtX ,double* pmtY,double* pmtZ,int location,float maxY,float maxZ)
{
	if(location==0)
	{
		//end caps are positioned at y position +/- half length of tank + radius
		//	cout<<" add 2* "<<maxY<<" to "<<*pmtY<<endl;
		*pmtY+=(maxY+maxZ);
	}
	if(location==2)
	{//		cout<<" subtract  2* "<<maxY<<" from "<<*pmtY<<endl;
		
		*pmtY=-(maxY+maxZ+*pmtY);
	}	
	if(location==1)
	{
		float angle=atan2(*pmtY,*pmtX)+(3.1415927/2.0);
		float rho=maxY*angle;
		*pmtX=rho;
		*pmtY=*pmtZ;		
	}
	*pmtZ=0.0;
	//cout<<" x,y,z of Unrolled phototube "<<*pmtX<<" "<<*pmtY<<" "<<*pmtZ<<endl;
	float xshift=maxY*(3.1415927);
	float x=*pmtX;
	if(x>xshift)x=x-(xshift*2);
	*pmtX=x;
	//cout<<" x,y,z of Unrolled phototube "<<*pmtX<<" "<<*pmtY<<" "<<*pmtZ<<endl;
	return;
}

TEveViewer* New2dView(TString name,TGLViewer::ECameraType type, TEveScene* scene)
{ 
	TEveViewer* View =gEve->SpawnNewViewer(name,name);
	View->AddScene(scene);  // add the special scene that only applies to this view 
	View->AddScene(gEve->GetGlobalScene()); // add the geometry information
	View->GetGLViewer()->SetCurrentCamera(type);
	View->GetGLViewer()->SetResetCamerasOnUpdate(kFALSE);
	return View;
}
void createPalettes()
{
	allPaletteHandlers.clear();
	pal3D = new TEveRGBAPalette();
	PaletteHandler* ph1 = new PaletteHandler("pal",pal3D); 
	pal3D->Connect("MinMaxValChanged()", "PaletteHandler", ph1,"PaletteChanged()");
	pal3D->SetLimits(0.0,1000.0);
	pal3D->SetMinMax(0.0,1000.0);
	pal3D->SetFixColorRange(kFALSE);
	pal3D->SetOverflowAction( TEveRGBAPalette::kLA_Clip);
	
	palUnrolled = new TEveRGBAPalette();//0, 4000.0);
	PaletteHandler* ph2 = new PaletteHandler("palUnrolled",palUnrolled); 
	palUnrolled->Connect("MinMaxValChanged()", "PaletteHandler", ph2,"PaletteChanged()");
	
	
	palUnrolled->SetLimits(0.0,1000.0);
	palUnrolled->SetMinMax(0.0,1000.0);
	palUnrolled->SetFixColorRange(kFALSE);
	palUnrolled->SetOverflowAction( TEveRGBAPalette::kLA_Clip);
}
//______________________________________________________________________________
void titus_eve_4Reco(TString inputFilename="",TString recoFileName="")
{
	//#include "titus_eve.h"
	
	// Main function, initializes the application.
	/*
	Open files
	*/
	const char *filetypes[] = {
		"ROOT files",    "*.root",    	
		"All files",     "*",
		0,               0
	};
	
	TString CurrentDirectory=gSystem->pwd();
	TString originalDirectory=CurrentDirectory;
	TGFileInfo fi;
	
	if(inputFilename.EqualTo(""))
	{
		fi.fFileTypes = filetypes;
		fi.fIniDir    = StrDup(CurrentDirectory);
		cout<<" Please choose your WCHsandbox  Titus file "<<endl;
		new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, &fi);
		if (!fi.fFilename) {
			cout<<" Nofile chosen "<<endl;
			return;
		}
		inputFilename=fi.fFilename;
	}
	cout<<" opening file "<<inputFilename<<endl;
	TFile* TitusFile = new TFile(inputFilename);
	
	TitusFile->ls();
	TitusEvents=(TTree *) TitusFile->Get("HitsTree");
	if(TitusEvents!=NULL)
		TitusPMTs=(TTree *) TitusFile->Get("PMTsTree");
	
	
	if(TitusEvents==NULL)
	{
		cout<<" No Titus Events TTree found in this file "<<endl;
		return;
	}
	TitusEvents->SetBranchAddress("evt", &evt, &b_evt);
	TitusEvents->SetBranchAddress("nhits", &nhits, &b_nhits);
	TitusEvents->SetBranchAddress("hit_time", hit_time, &b_hit_time);
	TitusEvents->SetBranchAddress("hit_x", hit_x, &b_hit_x);
	TitusEvents->SetBranchAddress("hit_y", hit_y, &b_hit_y);
	TitusEvents->SetBranchAddress("hit_z", hit_z, &b_hit_z);
	TitusEvents->SetBranchAddress("hit_PMTid", hit_PMTid, &b_hit_PMTid);
	TitusEvents->SetBranchAddress("npart", &npart, &b_npart);
	TitusEvents->SetBranchAddress("part_xStart", &part_xStart, &b_part_xStart);
	TitusEvents->SetBranchAddress("part_yStart", &part_yStart, &b_part_yStart);
	TitusEvents->SetBranchAddress("part_zStart", &part_zStart, &b_part_zStart);
	TitusEvents->SetBranchAddress("part_tStart", &part_tStart, &b_part_tStart);
	TitusEvents->SetBranchAddress("part_xEnd", &part_xEnd, &b_part_xEnd);
	TitusEvents->SetBranchAddress("part_yEnd", &part_yEnd, &b_part_yEnd);
	TitusEvents->SetBranchAddress("part_zEnd", &part_zEnd, &b_part_zEnd);
	TitusEvents->SetBranchAddress("part_tEnd", &part_tEnd, &b_part_tEnd);
	TitusEvents->SetBranchAddress("part_pxStart", &part_pxStart, &b_part_pxStart);
	TitusEvents->SetBranchAddress("part_pyStart", &part_pyStart, &b_part_pyStart);
	TitusEvents->SetBranchAddress("part_pzStart", &part_pzStart, &b_part_pzStart);
	TitusEvents->SetBranchAddress("part_pxEnd", &part_pxEnd, &b_part_pxEnd);
	TitusEvents->SetBranchAddress("part_pyEnd", &part_pyEnd, &b_part_pyEnd);
	TitusEvents->SetBranchAddress("part_pzEnd", &part_pzEnd, &b_part_pzEnd);
	TitusEvents->SetBranchAddress("part_KEstart", &part_KEstart, &b_part_KEstart);
	TitusEvents->SetBranchAddress("part_KEend", &part_KEend, &b_part_KEend);
	TitusEvents->SetBranchAddress("part_processStart", &part_processStart, &b_part_processStart);
	TitusEvents->SetBranchAddress("part_processEnd", &part_processEnd, &b_part_processEnd);
	TitusEvents->SetBranchAddress("part_parentid", &part_parentid, &b_part_parentid);
	TitusEvents->SetBranchAddress("part_trackid", &part_trackid, &b_part_trackid);
	TitusEvents->SetBranchAddress("part_pid", &part_pid, &b_part_pid);
	TitusEvents->SetBranchAddress("ncapturecount", &ncapturecount, &b_ncapturecount);
	TitusEvents->SetBranchAddress("neutroncount", &neutroncount, &b_neutroncount);
	TitusEvents->SetBranchAddress("capt_x", capt_x, &b_capt_x);
	TitusEvents->SetBranchAddress("capt_y", capt_y, &b_capt_y);
	TitusEvents->SetBranchAddress("capt_z", capt_z, &b_capt_z);
	TitusEvents->SetBranchAddress("capt_t0", capt_t0, &b_capt_t0);
	TitusEvents->SetBranchAddress("capt_E", capt_E, &b_capt_E);
	TitusEvents->SetBranchAddress("capt_num", capt_num, &b_capt_num);
	TitusEvents->SetBranchAddress("capt_pid", capt_pid, &b_capt_pid);
	TitusEvents->SetBranchAddress("capt_nucleus", capt_nucleus, &b_capt_nucleus);
	TitusEvents->SetBranchAddress("capt_nphot", capt_nphot, &b_capt_nphot);
	TitusEvents->SetBranchAddress("capt_ngamma", capt_ngamma, &b_capt_ngamma);
	TitusEvents->SetBranchAddress("mode", &mode, &b_mode);
	TitusEvents->SetBranchAddress("neutrino_E", &neutrino_E, &b_neutrino_E);
	TitusEvents->SetBranchAddress("neutrino_id", &neutrino_id, &b_neutrino_id);
	TitusEvents->SetBranchAddress("neutrino_px", &neutrino_px, &b_neutrino_px);
	TitusEvents->SetBranchAddress("neutrino_py", &neutrino_py, &b_neutrino_py);
	TitusEvents->SetBranchAddress("neutrino_pz", &neutrino_pz, &b_neutrino_pz);
	TitusEvents->SetBranchAddress("ntrks", &ntrks, &b_ntrks);
	TitusEvents->SetBranchAddress("nneutrons", &nneutrons, &b_nneutrons);
	TitusEvents->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
	TitusEvents->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
	TitusEvents->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
	TitusEvents->SetBranchAddress("mpid", &mpid, &b_mpid);
	TitusEvents->SetBranchAddress("px", &px, &b_px);
	TitusEvents->SetBranchAddress("py", &py, &b_py);
	TitusEvents->SetBranchAddress("pz", &pz, &b_pz);
	TitusEvents->SetBranchAddress("KE", &KE, &b_KE);
	
	if(TitusPMTs!=NULL)
	{
		TitusPMTs->SetBranchAddress("pmt_id", &pmt_id, &b_pmt_id);
		TitusPMTs->SetBranchAddress("pmt_pos_x", &pmt_pos_x, &b_pmt_pos_x);
		TitusPMTs->SetBranchAddress("pmt_pos_y", &pmt_pos_y, &b_pmt_pos_y);
		TitusPMTs->SetBranchAddress("pmt_pos_z", &pmt_pos_z, &b_pmt_pos_z);
		TitusPMTs->SetBranchAddress("pmt_size", &pmt_size, &b_pmt_size);
		TitusPMTs->SetBranchAddress("pmt_qe", &pmt_qe, &b_pmt_qe);
		TitusPMTs->SetBranchAddress("pmt_time_res", &pmt_time_res, &b_pmt_time_res);
		TitusPMTs->SetBranchAddress("is_lappd", &is_lappd, &b_is_lappd);
	}
	
	cout<<" Now open reconstructed file "<<endl;
	haveReconFile=kFALSE;
	if(recoFileName.EqualTo(""))
	{
		cout<<" Please choose the matching reconstructed file file "<<endl;
		new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, &fi);
		if (!fi.fFilename) {
			cout<<" Nofile chosen "<<endl;
		}
		else
		{
			recoFileName=fi.fFilename;
		}
	}
	if(!recoFileName.EqualTo(""))
	{
		cout<<" opening reconstruction file "<<recoFileName<<endl;
		TFile* TitusFile = new TFile(recoFileName);
		
		TitusFile->ls();
		Low_E=(TTree *) TitusFile->Get("Low_E");
		if(Low_E!=NULL)
			High_E_Electron=(TTree *) TitusFile->Get("High_E_Electron");
    if(High_E_Electron!=NULL)
			High_E_Muon=(TTree *) TitusFile->Get("High_E_Muon");
		if(High_E_Muon!=NULL)
			Final_Reconstruction=(TTree *) TitusFile->Get("Final_Reconstruction");
		
		if(Final_Reconstruction==NULL)
		{
			cout<<" No Final_Reconstruction  TTree found in this file "<<endl;
		}
		else
			haveReconFile=kTRUE;
	}
	
	if(Low_E!=NULL)
	{
		Low_E->SetBranchAddress("evt", &Low_E_evt, &b_Low_E_evt);
		Low_E->SetBranchAddress("Clusters", &Low_E_nClusters, &b_Low_E_nClusters);
		Low_E->SetBranchAddress("nSubevents", &Low_E_nSubevents, &b_Low_E_nSubevents);
		Low_E->SetBranchAddress("cluster", Low_E_cluster, &b_Low_E_cluster);
		Low_E->SetBranchAddress("recoVtxXLowE", recoVtxXLowE, &b_recoVtxXLowE);
		Low_E->SetBranchAddress("recoVtxYLowE", recoVtxYLowE, &b_recoVtxYLowE);
		Low_E->SetBranchAddress("recoVtxZLowE", recoVtxZLowE, &b_recoVtxZLowE);
		Low_E->SetBranchAddress("recoTimeLowE", recoTimeLowE, &b_recoTimeLowE);
		Low_E->SetBranchAddress("recoDirXLowE", recoDirXLowE, &b_recoDirXLowE);
		Low_E->SetBranchAddress("recoDirYLowE", recoDirYLowE, &b_recoDirYLowE);
		Low_E->SetBranchAddress("recoDirZLowE", recoDirZLowE, &b_recoDirZLowE);
		Low_E->SetBranchAddress("recoChkvAngleLowE", recoChkvAngleLowE, &b_recoChkvAngleLowE);
		Low_E->SetBranchAddress("recoEnergyLowE", recoEnergyLowE, &b_recoEnergyLowE);
	}
	if(High_E_Electron!=NULL)
	{
		
		High_E_Electron->SetBranchAddress("evt", &High_E_Electron_evt, &b_High_E_Electron_evt);
		High_E_Electron->SetBranchAddress("nClusters", &High_E_Electron_nClusters, &b_High_E_Electron_nClusters);
		High_E_Electron->SetBranchAddress("nSubevents", &High_E_Electron_nSubevents, &b_High_E_Electron_nSubevents);
		High_E_Electron->SetBranchAddress("cluster", High_E_Electron_cluster, &b_High_E_Electron_cluster);
		High_E_Electron->SetBranchAddress("ring", High_E_Electron_ring, &b_High_E_Electron_ring);
		High_E_Electron->SetBranchAddress("recoVtxXHighEElectron", recoVtxXHighEElectron, &b_recoVtxXHighEElectron);
		High_E_Electron->SetBranchAddress("recoVtxYHighEElectron", recoVtxYHighEElectron, &b_recoVtxYHighEElectron);
		High_E_Electron->SetBranchAddress("recoVtxZHighEElectron", recoVtxZHighEElectron, &b_recoVtxZHighEElectron);
		High_E_Electron->SetBranchAddress("recoTimeHighEElectron", recoTimeHighEElectron, &b_recoTimeHighEElectron);
		High_E_Electron->SetBranchAddress("recoDirXHighEElectron", recoDirXHighEElectron, &b_recoDirXHighEElectron);
		High_E_Electron->SetBranchAddress("recoDirYHighEElectron", recoDirYHighEElectron, &b_recoDirYHighEElectron);
		High_E_Electron->SetBranchAddress("recoDirZHighEElectron", recoDirZHighEElectron, &b_recoDirZHighEElectron);
		High_E_Electron->SetBranchAddress("recoChkvAngleHighEElectron", recoChkvAngleHighEElectron, &b_recoChkvAngleHighEElectron);
		High_E_Electron->SetBranchAddress("recoEnergyHighEElectron", recoEnergyHighEElectron, &b_recoEnergyHighEElectron);
		High_E_Electron->SetBranchAddress("recoLnLHighEElectron", recoLnLHighEElectron, &b_recoLnLHighEElectron);		
	}
	if(High_E_Muon!=NULL)
	{
		High_E_Muon->SetBranchAddress("evt", &High_E_Muon_evt, &b_High_E_Muon_evt);
		High_E_Muon->SetBranchAddress("nClusters", &High_E_Muon_nClusters, &b_High_E_Muon_nClusters);
		High_E_Muon->SetBranchAddress("nSubevents", &High_E_Muon_nSubevents, &b_High_E_Muon_nSubevents);
		High_E_Muon->SetBranchAddress("cluster", High_E_Muon_cluster, &b_High_E_Muon_cluster);
		High_E_Muon->SetBranchAddress("ring", High_E_Muon_ring, &b_High_E_Muon_ring);
		High_E_Muon->SetBranchAddress("recoVtxXHighEMuon", recoVtxXHighEMuon, &b_recoVtxXHighEMuon);
		High_E_Muon->SetBranchAddress("recoVtxYHighEMuon", recoVtxYHighEMuon, &b_recoVtxYHighEMuon);
		High_E_Muon->SetBranchAddress("recoVtxZHighEMuon", recoVtxZHighEMuon, &b_recoVtxZHighEMuon);
		High_E_Muon->SetBranchAddress("recoTimeHighEMuon", recoTimeHighEMuon, &b_recoTimeHighEMuon);
		High_E_Muon->SetBranchAddress("recoDirXHighEMuon", recoDirXHighEMuon, &b_recoDirXHighEMuon);
		High_E_Muon->SetBranchAddress("recoDirYHighEMuon", recoDirYHighEMuon, &b_recoDirYHighEMuon);
		High_E_Muon->SetBranchAddress("recoDirZHighEMuon", recoDirZHighEMuon, &b_recoDirZHighEMuon);
		High_E_Muon->SetBranchAddress("recoChkvAngleHighEMuon", recoChkvAngleHighEMuon, &b_recoChkvAngleHighEMuon);
		High_E_Muon->SetBranchAddress("recoEnergyHighEMuon", recoEnergyHighEMuon, &b_recoEnergyHighEMuon);
		High_E_Muon->SetBranchAddress("recoLnLHighEMuon", recoLnLHighEMuon, &b_recoLnLHighEMuon);
	}
	if(Final_Reconstruction!=NULL)
	{
		
		Final_Reconstruction->SetBranchAddress("evt", &Final_evt, &b_Final_evt);
		Final_Reconstruction->SetBranchAddress("nClusters", &Final_nClusters, &b_Final_nClusters);
		Final_Reconstruction->SetBranchAddress("nSubevents", &Final_nSubevents, &b_Final_nSubevents);
		Final_Reconstruction->SetBranchAddress("recoNRings", recoNRings, &b_recoNRings);
		Final_Reconstruction->SetBranchAddress("cluster", Final_cluster, &b_Final_cluster);
		Final_Reconstruction->SetBranchAddress("ring", Final_ring, &b_Final_ring);
		Final_Reconstruction->SetBranchAddress("recoVtxX", recoVtxX, &b_recoVtxX);
		Final_Reconstruction->SetBranchAddress("recoVtxY", recoVtxY, &b_recoVtxY);
		Final_Reconstruction->SetBranchAddress("recoVtxZ", recoVtxZ, &b_recoVtxZ);
		Final_Reconstruction->SetBranchAddress("recoTime", recoTime, &b_recoTime);
		Final_Reconstruction->SetBranchAddress("recoDirX", recoDirX, &b_recoDirX);
		Final_Reconstruction->SetBranchAddress("recoDirY", recoDirY, &b_recoDirY);
		Final_Reconstruction->SetBranchAddress("recoDirZ", recoDirZ, &b_recoDirZ);
		Final_Reconstruction->SetBranchAddress("recoChkvAngle", recoChkvAngle, &b_recoChkvAngle);
		Final_Reconstruction->SetBranchAddress("recoEnergy", recoEnergy, &b_recoEnergy);
		Final_Reconstruction->SetBranchAddress("recoPID", recoPID, &b_recoPID);
	}
	
	
	gSystem->cd(originalDirectory);
	
	
	/*
	Initialise Eve
	*/
	TEveManager::Create(kTRUE,"V");
	gEve->GetBrowser()->SetTabTitle("3D View",TRootBrowser::kRight,0);
	gEve->GetDefaultViewer()->SetName("3D View");
	gEve->GetDefaultViewer()->GetGLViewer()->SetResetCamerasOnUpdate(kFALSE);
	/*
	set up the geometry
	*/
	createGeometry(true);
	gGeoManager->SetTopVolume(SimpleVolume);
	gGeoManager->SetTopVisible();
	TGeoNode* node = gGeoManager->GetTopNode();
	TEveGeoTopNode* geomRoot = new TEveGeoTopNode(gGeoManager, node);
	geomRoot->SetVisLevel(4);
	geomRoot->GetNode()->GetVolume()->SetVisibility(kFALSE);
	gEve->AddGlobalElement(geomRoot); 
	/*
	Set up the event control GUI
	*/
	gEve->GetBrowser()->GetTabRight()->SetTab(1);
	make_gui();
	/*
	Initialise the unrolled geometry view
	*/
	gEve->GetBrowser()->GetTabRight()->SetTab(0);
	UnrolledScene = gEve->SpawnNewScene("Unrolled Event");
	flatGeometryScene = gEve->SpawnNewScene("Unrolled Geometry");
	UnrolledView = New2dView("Unrolled View",TGLViewer::kCameraOrthoXnOY,UnrolledScene);   	
	flatGeometryScene->AddElement(FlatGeometry);
	UnrolledView->AddScene(flatGeometryScene);
	TEveSceneInfo* gSI= (TEveSceneInfo*) (UnrolledView->FindChild("SI - Geometry scene"));
	if(gSI!=NULL)gSI->Delete();
	/*
	Initialise the palettes etc.
	*/
	gStyle->SetPalette(1, 0);
	
	
	TitusEvents->GetEntry(0) ;
	cout<<" There are "<<TitusEvents->GetEntries()<<" events "<<endl;
	/*
	load the first event
	*/
	load_event();
	/*
	Get Eve started
	*/
	gEve->GetDefaultGLViewer()->UpdateScene();
	gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.
	
}
/******************************************************************************/
// GUI
/******************************************************************************/
//______________________________________________________________________________
// 
// EvNavHandler class is needed to connect GUI signals.

class EvNavHandler
{
public:
	void Fwd()
	{
		if (event_id < (TitusEvents->GetEntries()-1)) {
			++event_id;
			load_event();
		} else {
			printf("Already at last event.\n");
		}
	}
	void Bck()
	{
		if (event_id > 0) {
			--event_id;
			load_event();
		} else {
			printf("Already at first event.\n");
		}
	}
	void GoToEventNumber(int)
	{
		Int_t eventRequested=NumberEntry->GetIntNumber();
		if(eventRequested>=TitusEvents->GetEntries())
		{
			cout<<" You  can not go beyond the  end of the file, there are "
			<<TitusEvents->GetEntries()<<" entries "<<endl;
			cout<<" numbered 0 to "<<TitusEvents->GetEntries()-1<<" ."<<endl;
			return;
		}
		if(eventRequested<0)
		{
			cout<<" Events numbers start at zero "<<endl;
			return;
		}
		event_id=eventRequested;
		load_event();
	}
	void SelecNEUTMode(int)
	{
		std::cout<<" you want to find an event with NEUT mode = "<<mode<<std::endl;
		Int_t mode=NEUTMode->GetIntNumber();
		std::cout<<" you want to find an event with NEUT mode = "<<mode<<std::endl;
		
	}
	void PrimaryOnly(Bool_t check)
	{
		drawPrimaryOnly=check;
		load_event();
	}
	
};


//______________________________________________________________________________
void make_gui()
{
	// Create minimal GUI for event navigation, etc.
	
	TEveBrowser* browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kLeft);
	
	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
	frmMain->SetWindowName("Titus Event Display");
	frmMain->SetCleanup(kDeepCleanup);
	TGCanvas*  fCanvasWindow = new TGCanvas(frmMain, 400, 240);
	TGCompositeFrame* fFrame = new TGCompositeFrame(fCanvasWindow->GetViewPort(), 10, 10, kVerticalFrame);
	fFrame->SetLayoutManager(new TGVerticalLayout(fFrame));    
	fCanvasWindow->SetContainer(fFrame);   
	// use hierarchical cleaning for container
	fFrame->SetCleanup(kDeepCleanup);  
	EvNavHandler* fh = new EvNavHandler;
	TGGroupFrame* Group;
	
	Group = new TGGroupFrame(fCanvasWindow->GetContainer(),"Event Navigation And Selection");
	TGHorizontalFrame* hf = new TGHorizontalFrame(Group,32,32);
	{
		
		TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
		TGPictureButton* b = 0;
		
		b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
		hf->AddFrame(b);
		b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");
		
		b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
		hf->AddFrame(b);
		b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
	}
	Group->AddFrame(hf);
	
	hf = new TGHorizontalFrame(Group,32,32);
	TGLayoutHints*  fButtonLayout =  new TGLayoutHints(kLHintsCenterY|kLHintsCenterX);
	{
		
		TGLabel* Label = new TGLabel(hf, "Event");
		hf->AddFrame(Label,fButtonLayout);
		NumberEntry = new TGNumberEntry(hf, 1, 7,4, TGNumberFormat::kNESInteger,
			TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMin,0);
		
		NumberEntry->Connect("ValueSet(Long_t)", "EvNavHandler", fh, "GoToEventNumber(int)");
		hf->AddFrame( NumberEntry, new TGLayoutHints(kLHintsCenterX, 5, 5, 5, 5));
		Label = new TGLabel(hf, " In File.");
		hf->AddFrame(Label,fButtonLayout);
	}
	Group->AddFrame(hf);
	hf = new TGHorizontalFrame(Group);
	{
		NEUTModeLabel = new TGLabel(hf, " NEUT mode: (unset) ");
		hf->AddFrame(NEUTModeLabel);
		//NEUTMode = new TGNumberEntry(hf, 1, 7,4, TGNumberFormat::kNESInteger,
		//	TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMin,0);
		//	NEUTMode->Connect("ValueSet(Long_t)", "EvNavHandler", fh, "SelecNEUTMode(int)");
		//hf->AddFrame( NumberEntry);//, new TGLayoutHints(kLHintsRightX, 5, 5, 5, 5));
	}
	Group->AddFrame(hf);
	hf = new TGHorizontalFrame(Group,32,32);
	{
		TGCheckButton *estat = new TGCheckButton(hf, "Primary MC tracks only",1);
		estat->SetToolTipText("Show just Primary MC tracks, or all");
		estat->SetState(kButtonDown);
		drawPrimaryOnly=kTRUE;
		hf->AddFrame(estat);
		estat->Connect("Toggled(Bool_t)", "EvNavHandler", fh, "PrimaryOnly(Bool_t)");
	}
	Group->AddFrame(hf);
	
	fCanvasWindow->AddFrame(Group);
	
	frmMain->AddFrame(fCanvasWindow,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
		0, 0, 2, 2));
	
	frmMain->MapSubwindows();
	frmMain->Resize();
	frmMain->MapWindow();
	
	browser->StopEmbedding();
	browser->SetTabTitle("Event Control", 0);
}
//______________________________________________________________________________
//void Init()
//{
//TitusEvents->SetBranchAddress("nhits", &nhits, &b_nhits);
//}
//______________________________________________________________________________

//______________________________________________________________________________


bool loadPMT(int hit_PMTid)
{
	for(int i=0;i<TitusPMTs->GetEntries();i++)
	{
		TitusPMTs->GetEntry(i);
		if(pmt_id==hit_PMTid)return true;
	}
	return false;
}
void load_event()
{
	createPalettes();
	
  
	// Load event specified in global event_id.
	// The contents of previous event are removed.
	
	//printf("Loading event %d.\n", event_id);
	
	gEve->GetViewers()->DeleteAnnotations();
	TEveEventManager* CurrentEvent =gEve->GetCurrentEvent();
	if(CurrentEvent != 0)CurrentEvent->DestroyElements();
	if( UnrolledScene !=0)UnrolledScene->DestroyElements();	
	
	TitusEvents->GetEvent(event_id);
	cout<<"Event : "<<evt;
	
	TEveBoxSet*  CherenkovHits= new TEveBoxSet("Hits ");
	CherenkovHits->SetPalette(pal3D);
	TEveBoxSet*  CherenkovHits2= new TEveBoxSet(Form("PMT Hits (Unrolled)"));
	CherenkovHits2->SetPalette(palUnrolled);
	
	
	
	CherenkovHits->Reset(TEveBoxSet::kBT_Cone, kFALSE, 64);
	CherenkovHits2->Reset(TEveBoxSet::kBT_Cone, kFALSE, 64);
	
	/*
	Draw the hit with colour at its correct location
	*/
	bool fDigitIsTime=kTRUE;
	float PMTRadius = pmt_size/2.0;
	int max=nhits;
	float minT=1e10;
	float maxT=-1e10;
	
	cout<<", NEUT  mode is:"<<mode<<endl;
	NEUTModeLabel->SetText(Form(" NEUT mode %i ",mode));
	std::map<int,hitStore> PMTmap;
	for(int i = 0;i<nhits;i++)
	{
		float Time=fmax(0.0,hit_time[i]); 
		minT = TMath::Min(minT,Time);
		maxT = TMath::Max(maxT,Time);
		// Add up all the hits contributing to a give PMT
		int tubeId=hit_PMTid[i];
		if(PMTmap.find(tubeId) != PMTmap.end())
		{
			PMTmap[tubeId].count++;
			PMTmap[tubeId].time+=Time;
		}
		else
		{
			hitStore temp;
			temp.count=1;
			temp.time=Time;
			PMTmap[tubeId]=temp;
		};
	}
	int maxCount=0;
	for(std::map<int,hitStore>::iterator pmt=PMTmap.begin();pmt!=PMTmap.end();pmt++)
	{
		if(pmt->second.count>maxCount)maxCount=pmt->second.count;
		pmt->second.time/=pmt->second.count;
	}
	/*
	Set range max to 90% of hit pmts
	*/
	int NPMTsHit=PMTmap.size();
	int limit=0.85*NPMTsHit;
	int lowLimit=0.02*NPMTsHit;
	
	std::set<float> orderedHits;
	for(std::map<int,hitStore>::iterator pmt=PMTmap.begin();pmt!=PMTmap.end();pmt++)
	{
		//	int key=pmt->first;
		hitStore value=pmt->second;
		//	pair<int,hitStore> p(key,value);
		orderedHits.insert(value.time);
	}
	//	for(std::set<pair<int, hitStore>, hitStoreCompare >::iterator pmt=orderedHits.begin();pmt!=orderedHits.end();pmt++)
  int count=0;
  float timeLimit=10000;
  bool foundLow=kFALSE;
  float timeMin=0.0;
  for(std::set<float>::iterator pmt=orderedHits.begin();pmt!=orderedHits.end();pmt++)
	{
		count++;
		if(!foundLow)
		{
			if(count>lowLimit)
			{
				timeMin=*pmt;
				foundLow=kTRUE;
			}
		}
		if(count>limit)
		{
			timeLimit=*pmt;
			break;
		}
	}
	timeLimit=fmin(timeLimit,150);
	if(fDigitIsTime)
	{
		if(maxT>100000)maxT=100000;
		pal3D->SetLimits(minT,maxT);
		pal3D->SetMin(timeMin);
		pal3D->SetMax(timeLimit);
		palUnrolled->SetLimits(minT,maxT);
		palUnrolled->SetMin(timeMin);
		palUnrolled->SetMax(timeLimit);
	}
	else
	{
		//pal->SetLimits(0,maxCharge);
		//pal->SetMin(0);
		//pal->SetMax(maxCharge);
	}	
	//for (int i = 0; i<max; i++){
	for(std::map<int,hitStore>::iterator pmt=PMTmap.begin();pmt!=PMTmap.end();pmt++)
	{
		//		int tubeId=hit_PMTid[i];
		int tubeId=pmt->first;
		//float Q=cDigiHit->GetQ(); 
		//float Time=hit_time[i]; 
		float Time=pmt->second.time; 
		int   nPhoton=pmt->second.count;
		if(!loadPMT(tubeId))continue;
		double pmtX = pmt_pos_x;
		double pmtY = pmt_pos_y;
		double pmtZ = pmt_pos_z;
		int location=PMTLocation();
		double pmtX2=pmtX;
		double pmtY2=pmtY;
		double pmtZ2=pmtZ;
		UnrollView(&pmtX2,&pmtY2,&pmtZ2,location,maxY,maxZ);
		
		
		float displayRadius=PMTRadius*(0.5*(1+(nPhoton/maxCount)));
		if(location==1)
		{
			float lengthXY=sqrt(pmtX*pmtX+pmtY*pmtY);
			float xd=pmtX/lengthXY;
			float yd=pmtY/lengthXY;
			CherenkovHits2->AddCone(TEveVector(pmtX2,pmtY2,pmtZ2),TEveVector(0.0,0.0,1.0) ,displayRadius ); 
			CherenkovHits->AddCone(TEveVector(pmtX,pmtY,pmtZ),TEveVector(xd,yd,0.0) ,displayRadius);
		}
		else
		{
			CherenkovHits2->AddCone(TEveVector(pmtX2,pmtY2,pmtZ2),TEveVector(0.0,0.0,1.0) ,displayRadius ); 
			CherenkovHits->AddCone(TEveVector(pmtX,pmtY,pmtZ),TEveVector(0.0,0.0,1.0) ,displayRadius ); 
		}
		if(fDigitIsTime)
		{
			CherenkovHits2->DigitValue(Time);
			CherenkovHits->DigitValue(Time);
		}
		else
		{
			//					CherenkovHits2->DigitValue(Q);
			//				CherenkovHits->DigitValue(Q);
		}
		
	}
	
	/*
	Plot the hit tubes 'hits'
	*/
	if(CherenkovHits->GetPlex()->VecSize()>0)
	{
		CherenkovHits->RefitPlex();
		CherenkovHits2->RefitPlex();
		TEveTrans& t = CherenkovHits2->RefMainTrans();
		t.SetPos(0.0,0.0,0.0);
		gEve->AddElement(CherenkovHits);
		UnrolledScene->AddElement(CherenkovHits2);
	}
	TEveElementList* Particles= new TEveElementList("Truth Particles");
	TEveElementList* NeutronCaptures= new TEveElementList("Neutron Captures");
	//cout<<" the neutrino has energy "<<neutrino_E<<" id "<<neutrino_id<<" x,y,z "<<neutrino_px
	//<<" "<<neutrino_py<<" "<<neutrino_pz<<" ntrks "<<ntrks<<" nneutrons "<<nneutrons<<" vtx "<<vtxx<<" "<<vtxy<<" "<<vtxz<<endl;
	
  // draw a line to represent the incoming neutrino
  TEveLine* neutrino = new TEveLine;
	neutrino->SetLineColor(kWhite);
  neutrino->SetLineWidth(2);
  neutrino->SetLineStyle(2);
  //cout            <<part_xEnd[i]<<" "<<part_yEnd[i]<<" "<<part_zEnd[i]<<endl;
	neutrino->SetNextPoint(10.*vtxx,10.*vtxy,10.*vtxz);
	// calculate start position at wall of cylinder
	// first calculate angle wrt z axis
	float theta = atan2(sqrt(neutrino_px*neutrino_px+neutrino_py*neutrino_py),neutrino_pz);
	float zToWall=(10.*vtxz)+maxZ;
	float length=zToWall;
	if(cos(theta)!=0)length=zToWall/cos(theta);
	neutrino->SetNextPoint(10.*vtxx-length*neutrino_px,10.*vtxy-length*neutrino_py,10.*vtxz-length*neutrino_pz);
	int pdgCode=neutrino_id;
	double ke=neutrino_E;
	TString Name;
	if(pdgCode==12)Name="electron anti-neutrino";
	if(pdgCode==-12)Name="electron anti-neutrino";
	if(pdgCode==14)Name=" muon neutrino";
	if(pdgCode==-14)Name="muon anti-neutrino";
	Name+=(Form(" energy %f MeV ",ke));
	
	if(abs(pdgCode)==12)neutrino->SetMainColor(kBlue);
	if(abs(pdgCode)==14)neutrino->SetMainColor(kGreen);
	neutrino->SetName(Name);
	neutrino->SetTitle(Name);
	Particles->AddElement(neutrino);
	
	// an array of pointers to TEveLine objects 
	TEveLine ** ls = new TEveLine*[npart];
	for(int i = 0;i<npart;i++)
	{
		ls[i] = new TEveLine(2); // create a new TEveLine for every true particle
		ls[i]->SetLineColor(kWhite);
		ls[i]->SetLineWidth(3);
		ls[i]->SetNextPoint(part_xStart[i],part_yStart[i],part_zStart[i]);
		ls[i]->SetNextPoint(part_xEnd[i],part_yEnd[i],part_zEnd[i]);
		int pdgCode=part_pid[i];
		if(pdgCode==0)continue;
		double ke=part_KEstart[i];
		TString Name(Form("MC particle %d, %f  MeV  ",i,ke));
		
		switch(pdgCode)
		{
		case 11: Name+=" (electron)"; break;
		case -11: Name+=" (positron)";break;
		case 12:Name+=" (electron neutrino)"; break;
		case -12 : Name+=" (electron anti-neutrino)";break;
		case 13 : Name+=" (muon)";break;
		case -13 :Name+="(anti-muon)";break;
		case 14:Name+=" (muon neutrino)"; break;
		case -14 :Name+=" (muon anti-neutrino)";break;
		case  2212 :Name+=" (proton)";break;
		case -2212 :Name+=" (anti-proton)";
		case 22 :Name+=" (photon)";break;
		case 211:Name+="(pi plus)";break;
		case -211:Name+="(pi minus)";break;
		case 111 :Name+="(pi0)";break;
		case 2112:Name+="(neutron)";break;
		default : Name+=Form(" pdg code %i",pdgCode);
		}
		
		if(abs(pdgCode)==11)ls[i]->SetMainColor(kYellow);
		if(abs(pdgCode)==12)ls[i]->SetMainColor(kBlue);
		if(abs(pdgCode)==13)ls[i]->SetMainColor(kMagenta);
		if(abs(pdgCode)==14)ls[i]->SetMainColor(kGreen);
		if(abs(pdgCode)==211)ls[i]->SetMainColor(kGreen);
		if(abs(pdgCode)==2212)ls[i]->SetMainColor(kRed);
		if(abs(pdgCode)==2112)ls[i]->SetMainColor(kBlue);
		ls[i]->SetName(Name);
		ls[i]->SetTitle(Name);
		
		if(part_parentid[i]==0) 
			Particles->AddElement(ls[i]); //primary particles are added to the top level 
		else
			if(drawPrimaryOnly)
			{
				/*
				If this is NOT a primary, and drawPrimaryOnly is TRUE, then set it to be not visible
				*/
				ls[i]->SetRnrSelf(kFALSE);
		  }
		  
		  
	}
	// now add other objects to their parent object to create a hierarchy in Eve
	std::map<int,TEveElement*> captureList;
	for(int i = 0; i<npart; i++){
		pdgCode=part_pid[i];
		if(part_parentid[i]==0) continue; // primaries are already there
		int parentIndex = 0;
		// look through the list of particles for the parent of this one
		// parentIndex is the pointer to this objects parent in the array
		while(part_trackid[parentIndex]!=part_parentid[i] && parentIndex<npart) parentIndex++;
		if(parentIndex==npart) Particles->AddElement(ls[i]); // hit the end without finding parent
		else 
		{
			if(pdgCode==2112)
			{
				// add the neutron to the list under its parent as usual
				ls[parentIndex]->AddElement(ls[i]); // add object to its parent
				// Look for neutron capture , add it to the neutron.
				TVector3 n(part_xEnd[i],part_yEnd[i],part_zEnd[i]);
				for(int iC = 0; iC<ncapturecount; iC++){
					TVector3 capture(capt_x[iC],capt_y[iC],capt_z[iC]);
					float diff=(capture-n).Mag();
					if(diff<0.0001)
					{
						TLorentzVector centre(capt_x[iC],capt_y[iC],capt_z[iC],capt_t0[iC]);
						float side=40;
						float vert[24];
						CreateCube(vert,centre,side);
						TEveBox* capture= new TEveBox(Form("Neutron Capture %i",iC));
						capture->SetVertices(vert);
						ls[i]->AddElement(capture); // apend a capture to the neutron
						captureList[i]=capture;	// remember this for later
						
						capture= new TEveBox(Form("Neutron Capture %i",iC));
						capture->SetVertices(vert);
						NeutronCaptures->AddElement(capture); // also maintain a separate list for convenience
						
					}
				}
			}
			else
				if(captureList.find(parentIndex)==captureList.end()) // if parent of this element is NOT a captured neutron
				ls[parentIndex]->AddElement(ls[i]); // add object to its parent
			else
				captureList[parentIndex]->AddElement(ls[i]); // add this element to the parent capture object
			
			
		}
	}
	
	gEve->AddElement(Particles);
	//	Now look at captures
	//	if(ncapturecount>0){
	/*		cout<<" there are "<<ncapturecount<<" neutron captures "<<endl;
	for(int i = 0; i<ncapturecount; i++){
	cout<<i<<" "<<capt_x[i]<<" "<<capt_y[i]<<" "<<capt_z[i]<<" "<<capt_t0[i]<<" "<<capt_num[i]
	<<" "<<capt_pid[i]<<" "<<capt_nucleus[i]<<" "<<capt_nphot[i]<<" "<<capt_ngamma[i]<<endl;
	// create a box around this point //
	TLorentzVector centre(capt_x[i],capt_y[i],capt_z[i],capt_t0[i]);
	float side=40;
	float vert[24];
	CreateCube(vert,centre,side);
	TEveBox* capture= new TEveBox(Form("Neutron Capture %i",i));
	capture->SetVertices(vert);
	NeutronCaptures->AddElement(capture);
	}
	*/
	if(ncapturecount>0)gEve->AddElement(NeutronCaptures);
	//}
	
  
	gEve->Redraw3D(kFALSE, kTRUE);
	NumberEntry->SetIntNumber(event_id);
	if(haveReconFile)loadReconstructedEvent();
}
void     CreateCube(float vert[24],TLorentzVector centreVec,float side){
	
	
	#define Xcoord 0
	#define Ycoord 1
	#define Zcoord 2
	
  
	float centre[3];
	centre[Xcoord]=centreVec.X();
	centre[Ycoord]=centreVec.Y();
	centre[Zcoord]=centreVec.Z();
	vert[Xcoord]=centre[Xcoord]-side/2;
	vert[Ycoord]=centre[Ycoord]-side/2;
	vert[Zcoord]=centre[Zcoord]-side/2;
	
	NextCubeVertex(&vert[0],side,Xcoord,1.0);
	NextCubeVertex(&vert[3],side,Ycoord,1.0);
	NextCubeVertex(&vert[6],side,Xcoord,-1.0);
	
	
	vert[12+Xcoord]=centre[Xcoord]-side/2.0;
	vert[12+Ycoord]=centre[Ycoord]-side/2.0;
	vert[12+Zcoord]=centre[Zcoord]+side/2.0;
	
	NextCubeVertex(&vert[12],side,Xcoord,1.0);
	NextCubeVertex(&vert[15],side,Ycoord,1.0);
	NextCubeVertex(&vert[18],side,Xcoord,-1.0);
	
	return ;
	
	
}
void               NextCubeVertex(float vert[6],float step,int changeMe,float sign){
	vert[3]=vert[0];
	vert[4]=vert[1];
	vert[5]=vert[2];
	vert[changeMe+3]+=sign*step;
	
}
void loadReconstructedEvent()
{
	//	cout<<" load reconstructed event here "<<endl;
	Final_Reconstruction->GetEvent(event_id);
	//	cout<<" Reconstructed Event : "<<Final_evt;
	//	cout<<" nclusters "<<Final_nClusters<<endl;
	//	for(int cluster=0;cluster<Final_nClusters;cluster++)
	//	{
	//		cout<<"cluster:"<<cluster<<" nrings "<<recoNRings[cluster]<<endl;
	//	}
	TEveElementList* ReconstructedObjects = new TEveElementList("Reconstructed Objects");
	for(int subevent=0;subevent<Final_nSubevents;subevent++)
	{
		
		// draw a line to represent the reconstructed particle
		if(!isnan(recoDirX[subevent]))
		{
			TEveLine* Track  = new TEveLine;
			Track->SetLineColor(kViolet);
			Track->SetLineWidth(4);
			Track->SetNextPoint(10.0*recoVtxX[subevent],10.0*recoVtxY[subevent],10.0*recoVtxZ[subevent]);
			float scale=1.0;
			float endX=10.0*recoVtxX[subevent]+scale*recoEnergy[subevent]*recoDirX[subevent];
			float endY=10.0*recoVtxY[subevent]+scale*recoEnergy[subevent]*recoDirY[subevent];
			float endZ=10.0*recoVtxZ[subevent]+scale*recoEnergy[subevent]*recoDirZ[subevent];
			Track->SetNextPoint(endX,endY,endZ);
			TString Name(Form("Reconstructed Track %i,  energy %f MeV ",subevent,recoEnergy[subevent]));
			Track->SetName(Name);
			Track->SetTitle(Name);
			ReconstructedObjects->AddElement(Track);
		}
		else
		{
			TLorentzVector centre(10.0*recoVtxX[subevent] ,10.0*recoVtxY[subevent] ,10.0*recoVtxZ[subevent], recoTime[subevent]);
			float side=40;
			float vert[24];
			CreateCube(vert,centre,side);
			TEveBox* RecoObject= new TEveBox(Form("Reconstructed Object  %i, energy: %f",subevent,recoEnergy[subevent]));
			RecoObject->SetVertices(vert);
			RecoObject->SetLineColor(kViolet);
			RecoObject->SetFillColor(kViolet);
			ReconstructedObjects->AddElement(RecoObject); // apend a capture to the neutron
		}
	}
	gEve->AddElement(ReconstructedObjects);
}
