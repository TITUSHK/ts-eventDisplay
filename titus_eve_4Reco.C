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
#include "TEveGeoShape.h"
#include "TEveGeoNode.h"

#include "TGLViewer.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include <map>


using namespace std;


const float tankR=5500;
const float tankZ=11000;

Double_t maxX,maxY,maxZ,minZ;

TTree * TitusEvents;
TTree * TitusPMTs;
void       make_gui();
void       load_event();

TGCheckButton *fColourIsTimeBox;
TGeoVolume *WorldVolume,*SimpleVolume;
TEveElementList *FlatGeometry;
TEveScene*  UnrolledScene;
TEveScene*  flatGeometryScene;
TEveViewer* UnrolledView;
TGNumberEntry* 		NumberEntry;
TGNumberEntry* 		NEUTMode;
TGLabel* NEUTModeLabel;
TEveViewer* New2dView(TString name,TGLViewer::ECameraType type, TEveScene* scene);
bool loadPMT(int hit_PMTid);

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

struct hitStore
{
	int count;
	float time;
}
;	
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
	cout<<" Phototube radius is "<<PMTRadius<<endl;
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

//______________________________________________________________________________
void titus_eve_4Reco()
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
	fi.fFileTypes = filetypes;
	fi.fIniDir    = StrDup(CurrentDirectory);
	cout<<" Please choose your WCHsandbox reconstructed Titus file "<<endl;
	new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, &fi);
	if (!fi.fFilename) {
		cout<<" Nofile chosen "<<endl;
		return;
	}
	cout<<" opening file "<<fi.fFilename<<endl;
	TFile* TitusFile = new TFile(fi.fFilename);
	
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
		/*	Int_t eventRequested=NumberEntry->GetIntNumber();
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
		*/
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
	EvNavHandler    *fh = new EvNavHandler;
	TGGroupFrame* Group;
	
	Group = new TGGroupFrame(fCanvasWindow->GetContainer(),"Event Navigation");
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
	
	
  
	// Load event specified in global event_id.
	// The contents of previous event are removed.
	
	//printf("Loading event %d.\n", event_id);
	
	gEve->GetViewers()->DeleteAnnotations();
	TEveEventManager* CurrentEvent =gEve->GetCurrentEvent();
	if(CurrentEvent != 0)CurrentEvent->DestroyElements();
	if( UnrolledScene !=0)UnrolledScene->DestroyElements();	
	
	TitusEvents->GetEvent(event_id);
	cout<<"Event : "<<evt<<endl;
	gStyle->SetPalette(1, 0);
	//TEveRGBAPalette* pal = new TEveRGBAPalette(150, 1000);
	TEveRGBAPalette* pal = new TEveRGBAPalette();//0, 4000.0);
	TEveRGBAPalette* palUnrolled = new TEveRGBAPalette();//0, 4000.0);
	
	TEveBoxSet*  CherenkovHits= new TEveBoxSet("Hits ");
	CherenkovHits->SetPalette(pal);
	TEveBoxSet*  CherenkovHits2= new TEveBoxSet(Form("PMT Hits (Unrolled)"));
	CherenkovHits2->SetPalette(palUnrolled);
	
	pal->SetLimits(0.0,500.0);
	pal->SetMinMax(0.0,500.0);
	pal->SetFixColorRange(kFALSE);
	pal->SetOverflowAction( TEveRGBAPalette::kLA_Clip);
	
	palUnrolled->SetLimits(0.0,500.0);
	palUnrolled->SetMinMax(0.0,500.0);
	palUnrolled->SetFixColorRange(kFALSE);
	palUnrolled->SetOverflowAction( TEveRGBAPalette::kLA_Clip);
	
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
	
	cout<<" mode :"<<mode<<endl;
	NEUTModeLabel->SetText(Form(" NEUT mode %i ",mode));
	std::map<int,hitStore> PMTmap;
	for(int i = 0;i<nhits;i++)
	{
		float Time=hit_time[i]; 
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
	if(fDigitIsTime)
	{
		if(maxT>1000)maxT=1000;
		pal->SetLimits(minT,maxT);
		pal->SetMin(minT);
		pal->SetMax(maxT);
		palUnrolled->SetLimits(minT,maxT);
		palUnrolled->SetMin(minT);
		palUnrolled->SetMax(maxT);
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
	// an array of pointers to TEveLine objects 
	TEveLine ** ls = new TEveLine*[npart];
	for(int i = 0;i<npart;i++)
	{
		ls[i] = new TEveLine(2); // create a new TEveLine for every true particle
		ls[i]->SetLineColor(kWhite);
		ls[i]->SetLineWidth(3);
		//cout<<" parton "<<part_xStart[i]<<" "<<part_yStart[i]<<" "<<part_zStart[i]<<" ";
		//cout            <<part_xEnd[i]<<" "<<part_yEnd[i]<<" "<<part_zEnd[i]<<endl;
		ls[i]->SetNextPoint(part_xStart[i],part_yStart[i],part_zStart[i]);
		ls[i]->SetNextPoint(part_xEnd[i],part_yEnd[i],part_zEnd[i]);
		//cout<<" pid "<<part_pid[i]<<endl;
		int pdgCode=part_pid[i];
		if(pdgCode==0)continue;
		double ke=part_KEstart[i];
		TString Name(Form("%i %f.*MeV PDG code %d ",(int)TMath::Max(0.,2-TMath::Log10(ke)),ke,pdgCode));
		
		
		if(pdgCode==11)Name+=" (electron)";
		if(pdgCode==-11)Name+=" (positron)";
		if(pdgCode==12)Name+=" (electron neutrino)";
		if(pdgCode==-12)Name+=" (electron anti-neutrino)";
		if(pdgCode==13)Name+=" (muon)";
		if(pdgCode==-13)Name+="(anti-muon)";
		if(pdgCode==14)Name+=" (muon neutrino)";
		if(pdgCode==-14)Name+=" (muon anti-neutrino)";
		if(pdgCode==2212)Name+=" (proton)";
		if(pdgCode==-2212)Name+=" (anti-proton)";
		if(pdgCode==22)Name+=" (photon)";
		
		if(pdgCode==211)Name+="pi plus";
		if(pdgCode==2112)Name+="neutron";
		
		if(abs(pdgCode)==11)ls[i]->SetMainColor(kYellow);
		if(abs(pdgCode)==12)ls[i]->SetMainColor(kBlue);
		if(abs(pdgCode)==13)ls[i]->SetMainColor(kMagenta);
		if(abs(pdgCode)==14)ls[i]->SetMainColor(kGreen);
		if(abs(pdgCode)==211)ls[i]->SetMainColor(kGreen);
		if(abs(pdgCode)==2212)ls[i]->SetMainColor(kRed);
		if(abs(pdgCode)==2112)ls[i]->SetMainColor(kBlue);
		ls[i]->SetName(Name);
		ls[i]->SetTitle(Name);
		if(part_parentid[i]==0) Particles->AddElement(ls[i]); //primary particles are added to the top level 
	}
	// now add other objects to their parent object to create a hierarchy in Eve
	for(int i = 0; i<npart; i++){
		if(part_parentid[i]==0) continue; // primaries are already there
		int parentIndex = 0;
		// look through the list of particles for the parent of this one
		// parentIndex is the pointer to this objects parent in the array
		while(part_trackid[parentIndex]!=part_parentid[i] && parentIndex<npart) parentIndex++;
		if(parentIndex==npart) Particles->AddElement(ls[i]); // hit the end without finding parent
		else ls[parentIndex]->AddElement(ls[i]); // add object to its parent
	}
	
	gEve->AddElement(Particles);
	//	TitusEvents->
	
	gEve->Redraw3D(kFALSE, kTRUE);
	NumberEntry->SetIntNumber(event_id);
}


