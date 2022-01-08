#include "ofApp.h"
#include "Fluid.h"


Fluid f(0.2, 0, 0.0000001);
int t = 0;
//--------------------------------------------------------------
void ofApp::setup(){
	ofSetBackgroundColor(0);
	ofSetFrameRate(60);
	ofSetWindowShape(N * SCALE, N * SCALE);
	ofSetWindowTitle("Fluid simulation");
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

	int cx = int(0.5 * ofGetWidth() / SCALE);
	int cy = int(0.5 * ofGetHeight() / SCALE);
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			f.add_density(cx + i, cy + j, ofRandom(100, 500));
		}
	}
	float angle = ofNoise(t) * TWO_PI;
	ofVec2f v = ofVec2f(1, 0).getRotated(angle);
	v *= (0.2);
	t += 0.01;
	f.add_velocity(cx, cy, 0, 0.1);

	f.step();
	f.render_d();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
