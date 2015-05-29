#!/bin/bash

ROOUNFOLD=RooUnfold-1.1.1
pushd $CMSSW_BASE/src/TTNJ/Unfold
if [ ! -d $ROOUNFOLD ]; then
  [ -f $ROOUNFOLD.tar.gz ] || wget http://hepunx.rl.ac.uk/~adye/software/unfold/$ROOUNFOLD.tar.gz
  tar xzf $ROOUNFOLD.tar.gz  
  pushd $ROOUNFOLD/src
  patch -R -p0 <<EOF
--- RooUnfoldTUnfold.cxx	2015-05-29 20:10:23.000000001 +0200
+++ RooUnfoldTUnfold.cxx.orig	2015-05-29 20:10:25.000000001 +0200
@@ -177,7 +177,7 @@
   // use automatic L-curve scan: start with taumin=taumax=0.0
   Double_t tauMin=0.0;
   Double_t tauMax=0.0;
-  //Int_t iBest;
+  Int_t iBest;
   TSpline *logTauX,*logTauY;
   TGraph *lCurve;
   // this method scans the parameter tau and finds the kink in the L curve
@@ -192,16 +192,14 @@
 #endif
   //_unf->SetConstraint(TUnfold::kEConstraintArea);
   if (!tau_set){
-    //iBest=_unf->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
-    _unf->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
+    iBest=_unf->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
     _tau=_unf->GetTau();  // save value, even if we don't use it unless tau_set
     cout <<"Lcurve scan chose tau= "<<_tau<<endl;
   }
   else{
     _unf->DoUnfold(_tau);
   }
-  TH1F* reco = new TH1F("h_rec", "reconstructed dist", _nt, 0, _nt);
-  _unf->GetOutput(reco);
+  TH1* reco=_unf->GetOutput("_rec","reconstructed dist",0,0);
   _rec.ResizeTo (_nt);
   for (int i=0;i<_nt;i++){
     _rec(i)=(reco->GetBinContent(i+1));
@@ -218,8 +216,7 @@
 {
   //Gets Covariance matrix
   if (!_unf) return;
-  TH2D* ematrix = new TH2D("ematrix", "error matrix", _nt, 0, _nt, _nt, 0, _nt);
-  _unf->GetEmatrix(ematrix);
+  TH2D* ematrix=_unf->GetEmatrix("ematrix","error matrix",0,0);
   _cov.ResizeTo (_nt,_nt);
   for (Int_t i= 0; i<_nt; i++) {
     for (Int_t j= 0; j<_nt; j++) {
EOF
  cd ..
  make -j8
  
  popd
fi

popd
