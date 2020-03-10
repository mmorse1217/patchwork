#include "face_map_obj.hpp"
#include "reflmap.h"
#include "mat2t.hpp"
#include <iomanip>

inline void Curvature( const Point3& F_u,  const Point3& F_v,
        const Point3& F_uu, const Point3& F_uv, const Point3& F_vv,
        const Point3& n,
        /*Point3& k1, Point3& k2,*/
        double& km1, double& km2 ) {
    double E,F,G,L,M,N;
    E = dot(F_u,F_u); F = dot(F_u,F_v); G = dot(F_v,F_v);
    L = dot(n,F_uu);  M = dot(n,F_uv);  N = dot(n,F_vv);
    // calculate the magnitude of the two curvatures
    // if the parameterization is nondegenerate,
    // F_u is not parallel to F_v, so EG = (|F_u||F_v|)^2 > (|F_u||F_v|cos(a))^2
    // = F^2 and EG - F^2 != 0
    double a = E*G-F*F, b = 2*F*M-E*N-L*G, c = L*N - M*M;
    // DZ make sure we  do not divide by zero
    //DZ!!  assert( a > DBL_EPSILON);
    km1 = (-b+sqrtf(max( b*b-4*a*c, 0.0)))/(2*a);
    km2 = (-b-sqrtf(max(b*b-4*a*c,0.0)))/(2*a);

    /*

       c = F*N-G*M; b = E*N-G*L; a = E*M-F*L;
       double x1,x2,y1,y2;
       if( fabs(a) > DBL_EPSILON) {
       y1 =  -b-sqrtf(max(b*b-4*a*c,0.0));  x1 = 2*a;
       y2 =  -b+sqrtf(max(b*b-4*a*c,0.0));  x2 = 2*a;
       } else if( fabs(c) > DBL_EPSILON ) {
    //DZ  if a = 0, c != 0,
    // replace at^2 + bt + c = 0 with a + bx + cx^2 = 0;
    // x and t are related as t_1 x_2 = 1,  t_2 x_2_1 = 1
    // compute sin and cos directly rather than tan
    x1 =  -b+sqrtf(max(b*b-4*a*c,0.0));  y1 = 2*c;
    x2 =  -b-sqrtf(max(b*b-4*a*c,0.0));  y2 = 2*c;
    } else if ( fabs(b) > DBL_EPSILON &&
    fabs(F) < DBL_EPSILON &&
    fabs(G) > DBL_EPSILON &&
    fabs(E) > DBL_EPSILON) {
    // a == c == 0 => either F != 0  && L = E*M/F && N = G*M/F
    // and it follows that dir. curvature is constant M/F
    // otherwise, F = 0 and M = 0
    // the curvature is (L*t^2 + N)/(E*t^2 + G)
    // with extremal values for t = 0 and t = infty
    // equal to N/G and L/E, E + G != 0 by assumption;
    // if G = 0, then N = 0 and dir. curvature is constant
    // if E = 0, then L = 0 and dir. curvature is constant
    // so the condition for the curvature to be nonconstant
    // is b != 0 && F == 0 && G != 0 && E != 0
    if( N/G > L/E) {
    km1 = N/G; y1 = 1.0; x1 = 0.0;
    km2 = L/E; y2 = 0.0; x2 = 1.0;
    } else {
    km2 = N/G; y2 = 1.0; x2 = 0.0;
    km1 = L/E; y1 = 0.0; x1 = 1.0;
    }
    } else { // constant dir. curvature
    // arbitrary values
    y1 = 1.0; x1 = 0.0;
    y2 = 0.0; x2 = 1.0;
    //    sin1 = 1.0; cos1 = 0.0;
    //    sin2 = 0.0; cos2 = 1.0;
    }
    Point3 Dk1 = F_u*y1 + F_v*x1;
    Point3 Dk2 = F_u*y2 + F_v*x2;
    // DZ normalize() is overly cautious; do normalization by hand
    // this is a hack, a better solution is needed
    double mag = Dk1.l2();
    if (mag != 0.0) Dk1 *= 1.0 / mag;
    mag = Dk2.l2();
    if (mag != 0.0) Dk2 *= 1.0 / mag;
    k1[0] = Dk1[0]; k1[1] = Dk1[1]; k1[2] = Dk1[2];
    k2[0] = Dk2[0]; k2[1] = Dk2[1]; k2[2] = Dk2[2];
    //DZ !!!  assert( km1 >= km2 );
    Point3 tmp = k1; k1 = k2; k2 = tmp;
    //  k1 = k1 * km1;
    //  k2 = k2 * km2;
    // use normalized curvatures

*/
}

//---------------------------------------------------
//  double zOffset = 1e-3;
void FaceMapSurfObj::renderObject()
{
    //---------- todo at first rendering
    int TTL = pow2(_lvl);
    int HLF = TTL/2;
    int NS = TTL ;	  // to render the half chart
    //do the following at the first rendering
    int chksz=8;
    Matrix<GLubyte> chkimage(chksz,chksz);
    if(_Vfxy.size()==0) {
        //1. make cubemap
        checkCubeMapExtension();  makeCubeMap(1024, 64);
        cubeMapOff();
        //      make1DTexMap(1024, 64); tex1DMapOff();
        //2. make checkboard
        for(int i=0; i<chksz; i++)
            for(int j=0; j<chksz; j++) {
                chkimage(i,j) = ((i%2==0)^(j%2==0)) * 255;
            }

        _texImage = NULL;
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glGenTextures(1, &_chkname);
        glBindTexture(GL_TEXTURE_2D, _chkname);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, chksz, chksz,	0, GL_LUMINANCE, GL_UNSIGNED_BYTE, chkimage.data());
        ifstream fin("COLORMAP"); if(fin.good()==false) { cerr<<"NO COLORMAP"<<endl; exit(0); }
        int nc; fin>>nc;
        _colormap.resize(nc,3);
        for(int i=0; i<nc; i++) {		
            fin>>_colormap(i,0);		
            fin>>_colormap(i,1);		
            fin>>_colormap(i,2);	 
        }
        fin.close();
        //3. generate data
        int numf = _bdsurf->numFaces();
        int numv = _bdsurf->numVertices();

        // MJM POSSIBLE BUG INTRODUCED HERE 
        // attempting to ignore outer for loop without rewriting all of the 
        // nested containers from hell...
        numv = 1;
        _Vfxy.resize(numv);	 _Vfpos.resize(numv);	 _Vfnor.resize(numv);
        _Vfd1.resize(numv);	 _Vfd2.resize(numv);	 _Vfd3.resize(numv);	 _Vfd4.resize(numv);	 _Vfd5.resize(numv);
        _Vfd6.resize(numv); 
        _Vfd7.resize(numv);
        _Vfd8.resize(numv);

        if (_gen == RENDER_EXTRA){//if gen ==2, render the whole chart
            NS = TTL-1;  
        }

        _Vfcgd.resize(numv);	 _Vfcps.resize(numv);
        _Vfcl0.resize(numv);	 _Vfcl1.resize(numv);

        if(_alt == EVAL_NORMAL || _alt == EVAL_HIGH_ORDER){
            cerr <<"generating data   "<< numv << ", " << numf << endl;	
            for(int V=0; V<numv; V++) {
                //int K = _bdsurf->valence(V);
                int K = numf;
                //if(K==0) continue;
                //if(K==4 && _gen == RENDER_EXTRA) continue; //only render extraordinary charts
                _Vfxy[V].resize(K);		_Vfpos[V].resize(K);		_Vfnor[V].resize(K);
                _Vfd1[V].resize(K);		_Vfd2[V].resize(K);		_Vfd3[V].resize(K);		  
                _Vfd4[V].resize(K);		_Vfd5[V].resize(K);
                _Vfd6[V].resize(K); _Vfd7[V].resize(K); _Vfd8[V].resize(K);

                for(int f=0; f<K; f++) {
                    Matrix<Point2>& cmxy  = _Vfxy[V][f];  cmxy.resize(NS+1,NS+1);
                    Matrix<Point3>& cmpos = _Vfpos[V][f]; cmpos.resize(NS+1,NS+1);
                    Matrix<Point3>& cmnor = _Vfnor[V][f]; cmnor.resize(NS+1,NS+1);
                    Matrix<double>& cmd1  = _Vfd1[V][f];  cmd1.resize(NS+1,NS+1);
                    Matrix<double>& cmd2  = _Vfd2[V][f];  cmd2.resize(NS+1,NS+1);
                    Matrix<double>& cmd3  = _Vfd3[V][f];  cmd3.resize(NS+1,NS+1);
                    Matrix<double>& cmd4  = _Vfd4[V][f];  cmd4.resize(NS+1,NS+1);
                    Matrix<double>& cmd5  = _Vfd5[V][f];  cmd5.resize(NS+1,NS+1);
                    Matrix<double>& cmd6  = _Vfd6[V][f];  cmd6.resize(NS+1,NS+1);
                    Matrix<double>& cmd7  = _Vfd7[V][f];  cmd7.resize(NS+1,NS+1);
                    Matrix<double>& cmd8  = _Vfd8[V][f];  cmd8.resize(NS+1,NS+1);
                    //cout << "blended face rendered: " << _bdsurf->_patches[f]._blended_face_to_sample<< endl;
                    double step = 1.0/double(TTL);
                    //cout << "num_samples_per_dim: " << NS << ", " << TTL << endl;
                    for(int j=0; j<=NS; j++)
                        for(int i=0; i<=NS; i++) {
                            //double cd[2];			
                            //cd[0] = i*step; cd[1] = j*step;
                            //double xy[2]; 
                            //_bdsurf->Vfcd2Vxy(BdSurf::EVAL_VALUE, V, f, cd, xy); 
                            //

                            double x_lower = _bdsurf->_patches[f]->_x_int.first;
                            double y_lower = _bdsurf->_patches[f]->_y_int.first;
                            double x_upper = _bdsurf->_patches[f]->_x_int.second;
                            double y_upper = _bdsurf->_patches[f]->_y_int.second;
                            //int ref_level = _bdsurf->_patches[f]->_level_of_refinement;
                            double xy[2];			
                            xy[0] = i*step; xy[1] = j*step;
                            if(i == 0 && j == 0){
                                xy[0] += 1e-7;
                                xy[1] += 1e-7;
                            }

                            Point3 ret[6];
                            
                            //Point3 blended_ret[6];
                            Point3* blended_ret;

                            if(       _gen== RENDER_FULL) {
                                _bdsurf->eval(  BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, f, xy, ret);
                                double cd[6];
                                //double st[6];

                                cd[0] =  xy[0]+1e-9; cd[1] = xy[1]+1e-9;
                                cd[0] = (x_upper - x_lower)*(cd[0]) + x_lower;
                                cd[1] = (y_upper - y_lower)*(cd[1]) + y_lower;
                                //if ( i == NS)
                                    //cd[0] -= 2e-14;
                                //if ( j == NS)
                                    //cd[1] -= 2e-14;
                                /*
                                cout << i << " , " << j << ": " << cd[0] << ", " << cd[1] << endl;

                                cout <<"x_int: ("<<  _bdsurf->_patches[f]._x_int.first << ", "
                                    << _bdsurf->_patches[f]._x_int.second <<")" <<   endl;

                                cout <<"y_int: ("<<  _bdsurf->_patches[f]._y_int.first << ", "
                                    << _bdsurf->_patches[f]._y_int.second <<")" <<   endl;
                                assert(cd[0] >_bdsurf->_patches[f]._x_int.first);
                                assert(cd[0] <=_bdsurf->_patches[f]._x_int.second);
                                assert(cd[1] >_bdsurf->_patches[f]._y_int.first);
                                assert(cd[1] <=_bdsurf->_patches[f]._y_int.second);
                                */
                                int blended_face_to_sample = (_bdsurf->_patches[f]->_blended_face_to_sample); // +3
                                //cout << blended_face_to_sample << ", " << f << endl;

                                //int local_vertex = 0;//(blended_face_to_sample)% 4; //1
                                //int global_vertex;
                                    //_bdsurf->_bdsurf->gpmesh().Fv2Vf(blended_face_to_sample, local_vertex).first;
                                //int local_face;
                                    //_bdsurf->_bdsurf->gpmesh().Fv2Vf(blended_face_to_sample, local_vertex).second;
                                    //
                                 /*
                                _bdsurf->_bdsurf->Fcd2Vfcd(
                                        BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV,
                                        blended_face_to_sample,
                                        cd,
                                        local_vertex,
                                        global_vertex,
                                        local_face,
                                        st
                                        );
                                for(int i =0 ; i < 6; i++)
                                    cd[i] = st[i];
                                _bdsurf->_bdsurf->Vfcd2Vxy(
                                        BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV,
                                        global_vertex,
                                        local_face,
                                        cd,
                                        st 
                                        );
                                _bdsurf->_bdsurf->eval(  BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, global_vertex, st, blended_ret);
                                */
                                //cout << blended_ret[0] << ", " << ret[0] << endl;
                                vector<Point3> position_and_derivs = 
                                    _bdsurf->evaluator()->evaluate(
                                            blended_face_to_sample, 
                                            BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, 
                                            Point2(cd));
                                blended_ret = position_and_derivs.data();
                                
                                //|BdSurf::EVAL_2ND_DERIV
                            } else if(_gen== RENDER_HALF || _gen == RENDER_EXTRA) {
                                _bdsurf->eval(  BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, f, xy, ret);
                                //_bdsurf->_bdsurf->eval(  BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, f, xy, blended_ret);
                                //_bdsurf->singleChartEval(BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, V, xy, ret);
                                //
                                //|BdSurf::EVAL_2ND_DERIV
                            }



                            //Point3 ret[6];
                            //-----------
                            cmxy(i,j) = Point2(xy[0], xy[1]);

                            cmpos(i,j) = ret[0];
                            cmnor(i,j) = cross(ret[1],ret[2]).dir();
                            cmd1(i,j) = sqrt(dot(ret[1],ret[1]) + dot(ret[2],ret[2]));
                            //cout << "fm: " << ret[0] << "blend: " <<  blended_ret[0] << endl;
                           //cout << "position error: " << (ret[0] - blended_ret[0]).l2() << endl;
                            cmd1(i,j) = (ret[0] - blended_ret[0]).l2()/blended_ret[0].l2();				  
                            Point3 n = cross(ret[1], ret[2]).dir();
                            Point3 bn = cross(blended_ret[1], blended_ret[2]).dir();
                            cmd2(i,j) = (n - bn).l2();///bn.l2();
                            				  cmd3(i,j) = 0;
                            cmd4(i,j) = 0;				  cmd5(i,j) = 0;
                            cmd6(i,j) = cmd7(i,j) = cmd8(i,j)= 0;
                            //cout << "face " << f<< " error:" << cmd1(i,j)<< endl;
                        }
                }
            } 


        }

        //_texture->initData(_Vfxy,_Vfpos,_Vfnor);
    }
    //  if(0){
    //  GpMesh& gpmesh = _bdsurf->gpmesh();
    double zOffset = 1e-3;
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_SURF) {
        //--------------------------
        glDepthRange(zOffset,1.0);//	 	 if(_mapOn==1) cubeMapOn();
        //glCullFace(GL_BACK);
        //glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        setMaterial();
        glColor3f(0.5,0.5,1);

        if(       _surfctrl==SURF_NONE) {
            //-----------------------------
            glEnable(GL_LIGHTING);
            for(int V=0; V<_Vfpos.size(); V++){
                for(int f=0; f<_Vfpos[V].size(); f++) {
                    Matrix<Point3>& cmpos = _Vfpos[V][f];
                    Matrix<Point3>& cmnor = _Vfnor[V][f];
                    for(int j=0; j<cmpos.n()-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<cmpos.m(); i++) {
                            glNormal3dv(cmnor(i,j+1));			  
                            glVertex3dv(cmpos(i,j+1));
                            glNormal3dv(cmnor(i,j));			
                            glVertex3dv(cmpos(i,j));
                        }
                        glEnd();
                    }
                }
            }
            glFlush();
        } else if(_surfctrl==SURF_CUBEMAP) {
            //-----------------------------
            glEnable(GL_LIGHTING);
            glMatrixMode(GL_TEXTURE);
            glScalef(0.5, 0.5, 1.0);
            cubeMapOn();
            //    tex1DMapOn();
            for(int V=0; V<_Vfpos.size(); V++)
                for(int f=0; f<_Vfpos[V].size(); f++) {
                    Matrix<Point3>& cmpos = _Vfpos[V][f];
                    Matrix<Point3>& cmnor = _Vfnor[V][f];
                    for(int j=0; j<cmpos.n()-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<cmpos.m(); i++) {
                            glNormal3dv(cmnor(i,j+1));				  glVertex3dv(cmpos(i,j+1));
                            glNormal3dv(cmnor(i,j));				  glVertex3dv(cmpos(i,j));
                        }
                        glEnd();
                    }
                }
            cubeMapOff();
            // tex1DMapOff();

            glLoadIdentity();
            glMatrixMode(GL_MODELVIEW);

            glFlush();
        } else if(_surfctrl==SURF_CHECKBOARD) {
            //-----------------------------
            glEnable(GL_LIGHTING);
            for(int V=0; V<_Vfpos.size(); V++)
                for(int f=0; f<_Vfpos[V].size(); f++) {
                    Matrix<Point3>& cmpos = _Vfpos[V][f];
                    Matrix<Point3>& cmnor = _Vfnor[V][f];
                    for(int j=0; j<HLF; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<HLF+1; i++) {
                            glNormal3dv(cmnor(i,j+1));				  glVertex3dv(cmpos(i,j+1));
                            glNormal3dv(cmnor(i,j));				  glVertex3dv(cmpos(i,j));
                        }
                        glEnd();
                    }
                }
            //render check board
            for(int i=0; i<chksz; i++)
                for(int j=0; j<chksz; j++) {
                    chkimage(i,j) = ((i%2==0)^(j%2==0)) * 255;
                }

            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glGenTextures(1, &_chkname);
            glBindTexture(GL_TEXTURE_2D, _chkname);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, chksz, chksz,	0, GL_LUMINANCE, GL_UNSIGNED_BYTE, chkimage.data());
            glEnable(GL_TEXTURE_2D);
            glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
            glDepthRange(zOffset/2, 1.0-zOffset/2);
            int stt = (_activevert==-1) ? 0 : _activevert;
            int end = (_activevert==-1) ? _Vfpos.size() : _activevert+1;
            for(int V=stt; V<end; V++)//int V = _activevert;
            for(int f=0; f<_Vfpos[V].size(); f++) {
                Matrix<Point2>& cmxy  = _Vfxy[V][f];
                Matrix<Point3>& cmpos = _Vfpos[V][f];
                Matrix<Point3>& cmnor = _Vfnor[V][f];
                for(int j=0; j<cmpos.n()-1; j++) {
                    glBegin(GL_QUAD_STRIP); //strip
                    for(int i=0; i<cmpos.m(); i++) {
                        glTexCoord2dv(cmxy(i,j+1));				  glNormal3dv(cmnor(i,j+1));				  glVertex3dv(cmpos(i,j+1));
                        glTexCoord2dv(cmxy(i,j));				  glNormal3dv(cmnor(i,j));				  glVertex3dv(cmpos(i,j));
                    }
                    glEnd();
                }
            }
            glDisable(GL_TEXTURE_2D);
            glDisable(GL_LIGHTING);
            /*
               glDepthRange(0,1.0-zOffset);
               glLineWidth(2.0);
               glColor3f(0.1,0.4,0.0);
               for(int f=0; f<_Vfpos[V].size(); f++) {
               Matrix<Point3>& cmpos = _Vfpos[V][f];
               glBegin(GL_LINE_STRIP);			 for(int i=0; i<cmpos.m(); i++)				glVertex3dv(cmpos(i,cmpos.n()-1));			 glEnd();
               glBegin(GL_LINE_STRIP);			 for(int j=0; j<cmpos.n(); j++)				glVertex3dv(cmpos(cmpos.m()-1,j));			 glEnd();
               }
               */
            glFlush();
        } else if(_surfctrl>=SURF_GD1 && _surfctrl<=SURF_GD8) {
            //-----------------------------
            glEnable(GL_LIGHTING);
            for(int V=0; V<_Vfpos.size(); V++)
                for(int f=0; f<_Vfpos[V].size(); f++) {
                    Matrix<Point3>& cmpos = _Vfpos[V][f];
                    Matrix<Point3>& cmnor = _Vfnor[V][f];
                    for(int j=0; j<HLF; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<HLF+1; i++) {
                            glNormal3dv(cmnor(i,j+1));				  glVertex3dv(cmpos(i,j+1));
                            glNormal3dv(cmnor(i,j));				  glVertex3dv(cmpos(i,j));
                        }
                        glEnd();
                    }
                }
            glDepthRange(zOffset/2, 1.0-zOffset/2);
            glDisable(GL_LIGHTING);
            vector< vector< Matrix<double> > >* dp = NULL;
            switch(_surfctrl) {
                case SURF_GD1: dp = &_Vfd1; cout<<"1 "; break;
                case SURF_GD2: dp = &_Vfd2; cout<<"2 "; break;
                case SURF_GD3: dp = &_Vfd3; cout<<"3 "; break;
                case SURF_GD4: dp = &_Vfd4; cout<<"4 "; break;
                case SURF_GD5: dp = &_Vfd5; cout<<"5 "; break;
                case SURF_GD6: dp = &_Vfd6; cout<<"6 "; break;
                case SURF_GD7: dp = &_Vfd7; cout<<"7 "; break;
                case SURF_GD8: dp = &_Vfd8; cout<<"8 "; break;
            }
            vector< vector< Matrix<double> > >& data = *dp;

            int stt = (_activevert==-1) ? 0 : _activevert;
            int end = (_activevert==-1) ? _Vfpos.size() : _activevert+1;
            //get min and max
            double mmin = SCL_MAX, mmax = -SCL_MAX;
            //int imax=-1, jmax=-1;
            Point3 maxpos;
            for(int V=stt; V<end; V++)
                for(int f=0; f<data[V].size(); f++) {			 //vector<NumMatrix> mags(data[V].size());//mags[f].resize(cm.m(), cm.n());
                    Matrix<Point3>& cmpos = _Vfpos[V][f];
                    Matrix<double>& cm = data[V][f]; //data
                    // cout<<cm.n()<<" "<<cm.m()<<endl;

                    for(int j=0; j<cm.n(); j++)
                        for(int i=0; i<cm.m(); i++) {
                            mmin = min(mmin, cm(i,j));

                            /* 
                               if(i!=0 && j!=0 && i!=cm.n()-1 && j!=cm.m()-1){

                               double test =( cm(i-1,j+1) + 2.*cm(i, j+1) + cm(i+1, j+1) +
                               2.*cm(i-1, j)  + 4.*cm(i,j) + 2.*cm(i+1, j)+
                               cm(i-1, j-1) + 2.*cm(i,j-1) + cm(i+1, j))*(1./16.);

                               mmax = max(mmax, test);
                               }
                               else*/

                            mmax = max(mmax, cm(i,j));
                            if(cm(i,j)==mmax) {maxpos = cmpos(i,j);}
                        }
                }

            //           cout<<setprecision(17);
            cout<<  mmin<< " " << mmax << endl; 
            //<< "  pos = "<<maxpos<<endl;		//mmin = 2*mmin - mmax; //make min = grey
            //      int nc = _colormap.m();
            for(int V=stt; V<end; V++)
                for(int f=0; f<data[V].size(); f++) {
                    Matrix<Point3>& cmpos = _Vfpos[V][f];			 //Matrix<Point3>& cmnor = _Vfnor[V][f];//Matrix<Point3>& cm = data[V][f];
                    Matrix<double>& cm = data[V][f]; //data
                    for(int j=0; j<cmpos.n()-1; j++) {			 //double tmp; Point3 clr;
                        glBegin(GL_QUAD_STRIP); //strip
                        for(int i=0; i<cmpos.m(); i++) {
                            {
                                double tmp = (cm(i,j+1)-mmin)/(mmax-mmin);
                                int idx = min(max((int)floor(tmp*64),0),63);
                                glColor3f(_colormap(idx,0), _colormap(idx,1), _colormap(idx,2));
                                glVertex3dv(cmpos(i,j+1));
                            }
                            {
                                double tmp = (cm(i,j  )-mmin)/(mmax-mmin);
                                int idx = min(max((int)floor(tmp*64),0),63);
                                glColor3f(_colormap(idx,0), _colormap(idx,1), _colormap(idx,2));
                                glVertex3dv(cmpos(i,j));
                            }
                        }
                        glEnd();
                    }
                }
            glDisable(GL_LIGHTING);
            glFlush();
        }
    }
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_FRAME) {
        //    cerr<<"FRAME"<<endl;
        glDepthRange(0,1.0-zOffset);
        glDisable(GL_LIGHTING);
        glColor3f(0.0, 0.0, 0.0);
        for(int V=0; V<_Vfpos.size(); V++) {
            for(int f=0; f<_Vfpos[V].size(); f++) {
                Matrix<Point3>& cmpos = _Vfpos[V][f];
                //	Matrix<Point3>& cmnor = _Vfnor[V][f];
                for(int j=0; j<cmpos.n(); j++) {
                    assert(!glCheck());
                    glBegin(GL_LINE_STRIP);
                    for(int i=0; i<cmpos.m(); i++) {
                        glVertex3dv(cmpos(i,j));
                    }
                    glEnd();
                }
                for(int i=0; i<cmpos.m(); i++) {
                    assert(!glCheck());
                    glBegin(GL_LINE_STRIP);
                    for(int j=0; j<cmpos.n(); j++) {
                        glVertex3dv(cmpos(i,j));
                    }
                    glEnd();
                }
            }
            glFlush();
        }
    }
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_CVTLINE) {
        cerr<<"CVTLINE"<<endl;
        glDepthRange(0,1.0-zOffset);
        glDisable(GL_LIGHTING);
        glColor3f(0,0,0);

        int NC = 32;

        int stt = (_activevert==-1) ? 0 : _activevert;
        int end = (_activevert==-1) ? _Vfpos.size() : _activevert+1;
        for(int V=stt; V<end; V++) {  //int V = _activevert;
            int K = _bdsurf->valence(V);
            if(K==0) continue;
            Matrix<bool>&   cmcgd = _Vfcgd[V];
            Matrix<Point3>& cmcps = _Vfcps[V];
            Matrix<Point3>& cmcl0 = _Vfcl0[V];
            Matrix<Point3>& cmcl1 = _Vfcl1[V];
            glBegin(GL_LINES);
            for(int j=0; j<=NC; j++)
                for(int i=0; i<=NC; i++)
                    if(cmcgd(i,j)==true) {
                        Point3 fm,to;
                        double h = 0.5;
                        fm = cmcps(i,j) - cmcl0(i,j)*h;				to = cmcps(i,j) + cmcl0(i,j)*h;				glVertex3dv(fm);				glVertex3dv(to);
                        fm = cmcps(i,j) - cmcl1(i,j)*h;				to = cmcps(i,j) + cmcl1(i,j)*h;				glVertex3dv(fm);				glVertex3dv(to);
                    }
            glEnd();
        }
        glFlush();
    }
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_INFLBDRY) {
        cerr<<"INFL"<<endl;
        glDepthRange(0,1.0-zOffset);
        glDisable(GL_LIGHTING);
        //    int numf = _bdsurf->numFaces();
        int numv = _bdsurf->numVertices();
        for(int V=0; V<numv; V++) {
            for(int f=0; f<_Vfpos[V].size(); f++) {
                Matrix<Point3>& cmpos = _Vfpos[V][f];
                //influence region boundary
                glColor3f(1.0, 0.0, 0.0);
                glBegin(GL_LINE_STRIP);			 for(int i=0; i<cmpos.m(); i++)				glVertex3dv(cmpos(i,cmpos.n()-1));			 glEnd();
                glBegin(GL_LINE_STRIP);			 for(int j=0; j<cmpos.n(); j++)				glVertex3dv(cmpos(cmpos.m()-1,j));			 glEnd();
            }
            glFlush();
        }
    }

    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_CPT) {
        glDisable(GL_LIGHTING);
        int nf = _bdsurf->gpmesh().numFaces();
        /*
        //face
        glDepthRange(zOffset,1.0);
        glColor3f(0.7, 0.7, 0.7);
        glBegin(GL_QUADS);

        for(int F=0; F<nf; F++) {
        for(int v=0; v<4; v++) {
        int V = _bdsurf->gpmesh().Fv2Vf(F,v).first;
        glVertex3dv(_bdsurf->gpmesh().vpoint(V).array());
        }
        }
        glEnd();
        */
        //line
        glDepthRange(0,1.0-zOffset);
        glColor3f(0.0, 0.0, 0.0);
        glLineWidth(2.0);
        for(int F=0; F<nf; F++) {
            glBegin(GL_LINE_LOOP);
            for(int v=0; v<4; v++) {
                int V = _bdsurf->gpmesh().Fv2Vf(F,v).first;
                glVertex3dv(_bdsurf->gpmesh().vpoint(V).array());
            }
            glEnd();
        }
        glFlush();
    }

    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_CCPTS) {
        glDisable(GL_LIGHTING);
        //int ctrllvl = _bdsurf->_bdsurf->subdivCtrlLevel();
        //ctrllvl = 5;

        //cout<<"level = "<<ctrllvl<<endl;
        glColor3f(1.0, 0.0, 0.0);
        glPointSize(3.0);

        //for(int F=0; F<_bdsurf->ccSurf().numFaces(); F++) { 
        for(int F=0; F<_bdsurf->numFaces(); F++) { 
            //CCRect<Point3>& cmpos =  _bdsurf->ccLimitPos()[F];
            // _bdsurf->ccSurf().pos(ctrllvl)[F];
            Vector<Point3> samples = _bdsurf->control_points(F);
            //vector<Point3> samples = (_bdsurf->samples(F));
                glBegin(GL_POINTS);
                //for(int i=0; i<samples.size(); i++) {
                for(int i=0; i<samples.length(); i++) {
                    glVertex3dv(samples(i).array());
                    glVertex3dv(samples(i).array());
                    //glVertex3dv(samples[i].array());
                    //glVertex3dv(samples[i].array());
                    //glVertex3dv(cmpos(i,j+1));
                    //glVertex3dv(cmpos(i,j));
                }
                glEnd();
            /*
            vector<Point3> samples = _bdsurf->samples(F);
                glBegin(GL_POINTS);
                //for(int i=0; i<samples.size(); i++) {
                for(int i=0; i<samples.size(); i++) {
                    glVertex3dv(samples[i].array());
                    glVertex3dv(samples[i].array());
                    //glVertex3dv(cmpos(i,j+1));
                    //glVertex3dv(cmpos(i,j));
                }
                glEnd();*/
            /*
            for(int j=0; j<pow2(ctrllvl); j++) {
                glBegin(GL_POINTS);
                for(int i=0; i<pow2(ctrllvl)+1; i++) {
                    glVertex3dv(cmpos(i,j+1));
                    glVertex3dv(cmpos(i,j));
                    //glVertex3dv(cmpos(i,j+1));
                    //glVertex3dv(cmpos(i,j));
                }
                glEnd();
            }*/

        }
        glFlush();
        setMaterial();
    }

    if(_renderctrl & RENDER_TEXTURE) {
        //    cerr<<"TEXTURE"<<endl;
        glEnable(GL_LIGHTING);
        //-----------------------------
        if(_texImage == NULL){
            glClearColor( 0.0, 0.0, 0.0, 0.0);
            _texImage = _texture->read(&_texSize[0], &_texSize[1]);
            cerr << _texSize[0] << "\t" << _texSize[1] << "\n";
            gluBuild2DMipmaps(GL_TEXTURE_2D, 4, _texSize[0], _texSize[1], GL_RGBA, GL_UNSIGNED_BYTE, _texImage);
            textureRenderType = GL_CLAMP;

            _texture->calcDisplacementMap();
            const GLfloat border_color[4] = {0.,0.,0.,0.};
            glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, &border_color[0]);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureRenderType);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureRenderType);
            glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
        }
        /*
         * render texture
         */
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureRenderType);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureRenderType);

        glEnable(GL_TEXTURE_2D);
        glDepthRange(zOffset/2, 1.0-zOffset/2);

        /*
         * calculate texture coordinates
         */
        if(updateTexture){
            _texture->updateTextureParam();
            updateTexture = false;
        }

        int nfaces = _texture->numFaces();
        for(int ii=0; ii<nfaces; ii++){
            Matrix<Point2>& cmxy  = _texture->Vfxy(ii);
            Matrix<Point3>& cmpos = _texture->Vfpos(ii);
            Matrix<Point3>& cmnor = _texture->Vfnor(ii);

            int mmax = cmpos.m();
            int nmax = cmpos.n();

            if(_renderctrl & RENDER_SURF){
                for(int j=0; j<nmax-1; j++) {

                    glBegin(GL_QUAD_STRIP); //strip

                    for(int i=0; i<mmax; i++) {
                        glTexCoord2dv(cmxy (i,j+1));
                        glNormal3dv  (cmnor(i,j+1));
                        glVertex3dv  (cmpos(i,j+1));
                        glTexCoord2dv(cmxy (i,j  ));
                        glNormal3dv  (cmnor(i,j  ));
                        glVertex3dv  (cmpos(i,j  ));
                    }
                    glEnd();
                }
            }
            if(_renderctrl & RENDER_FRAME) {
                //    cerr<<"FRAME"<<endl;
                glDepthRange(0,1.0-zOffset);
                glDisable(GL_LIGHTING);
                glColor3f(0.0, 0.0, 0.0);
                for(int j=0; j<nmax; j++) {
                    assert(!glCheck());
                    glBegin(GL_LINE_STRIP);
                    for(int i=0; i<mmax; i++) {
                        glVertex3dv(cmpos(i,j));
                    }
                    glEnd();
                }
                for(int i=0; i<mmax; i++) {
                    assert(!glCheck());
                    glBegin(GL_LINE_STRIP);
                    for(int j=0; j<nmax; j++) {
                        glVertex3dv(cmpos(i,j));
                    }
                    glEnd();
                }
                glEnable(GL_LIGHTING);
                glDepthRange(zOffset/2, 1.0-zOffset/2);
            }
        }
        glDisable(GL_TEXTURE_2D);    glDisable(GL_LIGHTING);
        glFlush();
    }
}


void FaceMapSurfObj::render(){
    if(selectionMode){
        doSelection();
    }
    renderObject();
}

void FaceMapSurfObj::doSelection(){
    int f;
    Point2 c;
    Point3 sc(0,0,0);
    Point3 sn(0,0,0);

    if(_movementctrl == TEX_MOVE_COLOR){
        _texture->storeOldPos();
        if(selectColor(mousePosition[0], mousePosition[1], f, c, sc, sn)){
            _texture->storePos(f,c,sc,sn);
            updateTexture = true;
            selectionMode = false;
        }
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    }else{
        selectSurfacePoint(mousePosition[0] , mousePosition[1], /*f, c,*/ sc, sn);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    }
}

//---------------------------------------------------
void FaceMapSurfObj::key(unsigned char k)
{
    //LEXING: q,c,w are reserved
    switch(k) {
        case 'r':
            _renderctrl = _renderctrl ^ RENDER_SURF; 	 break;
        case 'f':
            _renderctrl = _renderctrl ^ RENDER_FRAME; 	 break;
        case 'i':
            _renderctrl = _renderctrl ^ RENDER_INFLBDRY;	 break;
        case 'p':
            _renderctrl = _renderctrl ^ RENDER_CPT;        break;
        case 'k':
            _renderctrl = _renderctrl ^ RENDER_CCPTS;        break;
        case 'l':
            _renderctrl = _renderctrl ^ RENDER_CVTLINE;	 break;
        case 's':
            _surfctrl = (_surfctrl+1) % SURF_TTL; break;
        case 'v':
            if(     _activevert==-1) _activevert=0;
            else if(_activevert==_bdsurf->numVertices()-1) _activevert=-1;
            else _activevert++;
            //cout<<"active_vert = "<<_activevert<<endl;
            break;  //_activevert = (_activevert+1) % _bdsurf->numVertices(); break;
        case 'h':
            printf("r :Render Surface On| Off\n");
            printf("f: Show Wireframe On|Off\n");
            printf("i: Show Influence Boundary On|Off\n");
            printf("p: Show Control Points On|Off\n");
            printf("s: Loop in Surface Mode\n");
            printf("t: Loop Texture Mode texture|displacemen map\n");
            printf("a: Rotate texture CCW\n");
            printf("d: Rotate texture CW\n");
            printf("z: Switch to GL_CLAMP\n");
            printf("x: Switch to GL_CLAMP_TO_BORDER\n");
            printf("+: Increase texture size\n");
            printf("-: Decrease texture size\n");
            printf("]: Increase displacement map size\n");
            printf("[: Decrease displacement map size\n");
            break;
        case 't':
            _texture->switchTextureMode();
            _texImage = NULL;
            updateTexture = true;
            break;
        case 'a':
            _texture->rotateTexture(ROTATE_CCW);
            updateTexture = true;
            break;
        case 'd':
            _texture->rotateTexture(ROTATE_CW);
            updateTexture = true;
            break;
        case 'z':
            textureRenderType = GL_CLAMP;
            break;
        case 'x':
            textureRenderType = GL_CLAMP_TO_BORDER;
            break;
        case '-':
            _texture->scaleTexture(-0.05);
            updateTexture = true;
            break;
        case '+':
        case '=':
            _texture->scaleTexture(0.05);
            updateTexture = true;
            break;
        case '[':
            _texture->scaleDispMap(-0.05);
            updateTexture = true;
            break;
        case ']':
            _texture->scaleDispMap(0.05);
            updateTexture = true;
            break;
    }
    glutPostRedisplay();
}

void FaceMapSurfObj::specialKey(unsigned char k){
    switch(k){
        case GLUT_KEY_LEFT:
            _texture->moveTexture(MOVE_LEFT);
            updateTexture = true;
            break;
        case GLUT_KEY_RIGHT:
            _texture->moveTexture(MOVE_RIGHT);
            updateTexture = true;
            break;
        case GLUT_KEY_UP:
            _texture->moveTexture(MOVE_UP);
            updateTexture = true;
            break;
        case GLUT_KEY_DOWN:
            _texture->moveTexture(MOVE_DOWN);
            updateTexture = true;
            break;
        default:
            break;
    }
}


void FaceMapSurfObj::mouse(int button, int state, int x, int y){
    //  cout << x <<"\t" << y << endl;
    if(state == GLUT_DOWN){
        if(button == GLUT_LEFT_BUTTON){
            mouseOldPosition[0] = x;
            mouseOldPosition[1] = y;
        }
    }
    if(state == GLUT_DOWN && button == GLUT_LEFT_BUTTON){
        int f = 0;
        Point2 c;
        Point3 sc(0,0,0);
        Point3 sn(0,0,0);
        int modKey = glutGetModifiers();
        if (modKey == GLUT_ACTIVE_CTRL){
            _renderctrl = _renderctrl ^ RENDER_TEXTURE;
            if(_renderctrl & RENDER_TEXTURE){
                selectSurfacePoint(mouseOldPosition[0], mouseOldPosition[1], 
                        /*f, c,*/ sc, sn);
                _texture->setTextureMode(f, c, sc, sn);
                int power = pow2(_lvl);
                int n = max(_bdsurf->numFaces(),power);
                frexp((double)n, &power);
                base = pow2(power);
                updateTexture = true;
                selectionMode = false;
            }else{
                _texture->remove();
                _texImage = NULL;
            }
        }
    }
}

void FaceMapSurfObj::motion(int x, int y){

    mousePosition[0] = x;
    mousePosition[1] = y;

    double fw = glutGet(GLUT_SCREEN_WIDTH)/4.;
    double d[2] = {(x-mouseOldPosition[0])/fw, (y-mouseOldPosition[1])/fw};

    if(_renderctrl & RENDER_TEXTURE){
        switch(_movementctrl){
            case TEX_MOVE_UV:
                _texture->moveTextureCenter(d);
                updateTexture = true;
                break;
            case TEX_MOVE_SELECT:
                if(abs(x - mouseOldPosition[0]) > 1 || abs(y - mouseOldPosition[1]) > 1){
                    selectionMode = true;
                }
                break;
            case TEX_MOVE_COLOR:
                selectionMode = true;
                break;
        }
    }
    mouseOldPosition[0] = x;
    mouseOldPosition[1] = y;
}

//---------------------------------------------------
void FaceMapSurfObj::setMaterial()
{
    if(0) {
        float matSpecular[] = {0.478814, 0.457627, 0.5, 1.0};
        float matDiffuse[] =  {0.15, 0.452647, 0.154303, 1,0};
        float matAmbient[] = {0.15, 0.452647, 0.154303, 1.0};
        float matShininess = 10.0;
        float matBackSpecular[] = {0.5, 0.1, 0.75, 1.0};
        float matBackAmbient[] =  {0.15, 0.154303, 0.754303, 1.0};
        float matBackDiffuse[] =  {0.15, 0.154303, 0.754303, 1.0};
        float matBackShininess = 10.0;
        glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
        glMaterialfv(GL_FRONT,  GL_AMBIENT, matAmbient);
        glMaterialf(GL_FRONT, GL_SHININESS, matShininess);
        glMaterialfv(GL_BACK, GL_DIFFUSE, matBackDiffuse);
        glMaterialfv(GL_BACK, GL_SPECULAR, matBackSpecular);
        glMaterialfv(GL_BACK, GL_AMBIENT, matBackAmbient);
        glMaterialf(GL_BACK, GL_SHININESS, matBackShininess);
    }
    else{
        //    if(0) { //HENNING's material
        float matSpecular[] = {0.478814, 0.457627, 0.5};

        float matAmbient[] =  {0.75, 0.652647, 0.154303};
        float matDiffuse[] =  {0.75, 0.652647, 0.154303};



        float matFrontShininess = 25.0;

        float matBackSpecular[] = {0.1596,   0.1525,   0.1667};
        float matBackAmbient[] =  {0.3750,   0.3263,   0.0772};
        float matBackDiffuse[] =  {0.3750,   0.3263,   0.0772};
        float matBackShininess = 100.0;

        glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
        glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbient);
        glMaterialf(GL_FRONT, GL_SHININESS, matFrontShininess);


        /*      glMaterialfv(GL_FRONT, GL_DIFFUSE, matBackDiffuse);
                glMaterialfv(GL_FRONT, GL_SPECULAR, matBackSpecular);
                glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbient);
                glMaterialf(GL_FRONT, GL_SHININESS, matFrontShininess);
                */
        glMaterialfv(GL_BACK, GL_DIFFUSE, matBackDiffuse);
        glMaterialfv(GL_BACK, GL_SPECULAR, matBackSpecular);
        glMaterialfv(GL_BACK, GL_AMBIENT, matBackAmbient);
        glMaterialf(GL_BACK, GL_SHININESS, matBackShininess);	 
    }  
    glCheck();
    }


#define BUFSIZE 51200
    void FaceMapSurfObj::selectSurfacePoint( int x, int y, /*int& face, double cd[2],*/ Point3& p, Point3& n){

        glCullFace(GL_BACK);
        glEnable(GL_CULL_FACE);

        /*
         * make selection using depth buffer
         */
        GLuint selectBuf[BUFSIZE];
        GLint hits;
        GLint viewport[4];

        drawRects(GL_RENDER);
        glGetIntegerv(GL_VIEWPORT, viewport);

        glSelectBuffer(BUFSIZE, selectBuf);
        (void) glRenderMode(GL_SELECT);

        glInitNames();
        glPushName(0);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        /*  create 5x5 pixel picking region near cursor location */
        //cout << x << "\t" << y << "\t" << viewport[3] << endl;
        gluPickMatrix((GLdouble) x, (GLdouble) (y),
                1.0, 1.0, viewport);

        gluPerspective(fovy, aspect, znear, zfar);
        drawRects(GL_SELECT);
        glPopMatrix();

        hits = glRenderMode(GL_RENDER);
        GLuint quadNo = processHits(hits, selectBuf);

        // get face number and location
        int TTL = pow2(_lvl);
        int HLF = TTL/2;
        int count = 1;
        //int Fnum = 0;
        int Vnum = 0;
        int fnum = 0;
        int iindex = 0;
        int jindex = 0;
        if(hits!= 0){
            int nfaces = _texture->numFaces();
            if(nfaces == 0){
                for(int V=0; V<_Vfpos.size(); V++){
                    for(int f=0; f<_Vfpos[V].size(); f++) {
                        for(int j=0; j<HLF; j++) {
                            for(int i=0; i<HLF; i++) {
                                if(count == quadNo){
                                    Vnum = V;
                                    fnum = f;
                                    iindex = i;
                                    jindex = j;
                                }
                                count ++;
                            }
                        }
                    }
                }
            }else{
                //Fnum = (quadNo >> (power*2)) & (base-1);
                iindex = (quadNo >> power) & (base-1);
                jindex = (quadNo) & (base-1);
            }

            if(iindex == 0 && jindex ==0){
                iindex++;
            }

            //int nsteps = TTL;
            //double step = 1 / double(nsteps);        //defined in a unit square, length of every interval
            //double lcd[2];
            //lcd[0] = iindex*step;
            //lcd[1] = jindex*step;

            if(nfaces ==0){
                //_bdsurf->Vfcd2Fcd(BdSurf::EVAL_VALUE, Vnum, fnum, lcd, face, cd);
                //cd = lcd;

                p = _Vfpos[Vnum][fnum](iindex, jindex);

                n = _Vfnor[Vnum][fnum](iindex, jindex);
                //        }else{
                //            face = _textureFaces[Fnum];
                //            cd[0] = lcd[0];
                //            cd[1] = lcd[1];
                //            p = _Vftpos[Fnum](iindex,jindex);
                //            n = _Vftnor[Fnum](iindex,jindex);
        }
        }
        glDisable(GL_CULL_FACE);
    }

    bool FaceMapSurfObj::selectColor(int x, int y, int& face, double cd[2], Point3& p, Point3& n){

        int Fnum = 0;
        int iindex = 0;
        int jindex = 0;

        pickRender();

        GLubyte buff[10];

        glReadPixels(x,y,1,1,GL_RGB,GL_UNSIGNED_BYTE, buff);
        int r = buff[0]; int g = buff[1]; int b = buff[2];
        //    cout << " r : " << r << " g : " << g << " b : " << b << " \n";
        //    cout << " face : " << face<< "\n";
        if( r == 255 && g == 255 && b == 255){
            //        cout << " face : " << face<< "\n";
            return false;
        }
        int num = r*256*256 + g*256 + b;
        vector<Matrix<Point3> >& Vpos = _texture->Vftpos();
        vector<Matrix<Point3> >& Vnor = _texture->Vftnor();
        int nTexFaces = _texture->numFaces();
        int count =0;
        for(int ii=0; ii<nTexFaces; ii++){
            Matrix<Point3>& cmpos = Vpos[ii];
            int mmax = cmpos.m();
            int nmax = cmpos.n();
            for(int j=0; j<nmax-1; j++) {
                for(int i=0; i<mmax-1; i++) {
                    if(count == num){
                        Fnum = ii;
                        iindex = i;
                        jindex = j;
                    }
                    count++;
                }
            }
        }
        if(iindex == 0 && jindex ==0){
            iindex++;
        }
        int nsteps = pow2(_lvl);
        double step = 1 / double(nsteps);        //defined in a unit square, length of every interval
        double lcd[2];
        lcd[0] = iindex*step;
        lcd[1] = jindex*step;
        face = _texture->textFace(Fnum);
        cd[0] = lcd[0];
        cd[1] = lcd[1];
        p = Vpos[Fnum](iindex,jindex);
        n = Vnor[Fnum](iindex,jindex);
        glFlush();

        return true;
    }

    void FaceMapSurfObj::pickRender(){

        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        glDisable(GL_DITHER);
        glDisable(GL_LIGHTING);

        int count = 0;
        int rgb[3];
        int nfaces = _texture->numFaces();
        for(int ii=0; ii<nfaces; ii++){
            Matrix<Point3>& cmpos = _texture->Vftpos(ii);
            int mmax = cmpos.m();
            int nmax = cmpos.n();
            glBegin(GL_QUADS);
            for(int j=0; j<nmax-1; j++) {
                for(int i=0; i<mmax-1; i++) {
                    unpack(count, rgb);
                    glColor3ub((GLubyte)rgb[0], (GLubyte)rgb[1], (GLubyte)rgb[2]);
                    glVertex3dv(cmpos(i,j+1));
                    glVertex3dv(cmpos(i,j));
                    glVertex3dv(cmpos(i+1,j));
                    glVertex3dv(cmpos(i+1,j+1));
                    count++;
                }
            }
            glEnd();
        }
    }

    void  FaceMapSurfObj::unpack(int v, int r[3]){
        r[0] = (v >> 16) & 255;
        r[1] = (v >> 8) & 255;
        r[2] = (v) & 255;
    }

    void FaceMapSurfObj::drawRects(GLenum mode){
        int TTL = pow2(_lvl);
        int HLF = TTL/2;
        GLuint count = 1;

        int nfaces = _texture->numFaces();
        if(nfaces == 0){
            for(int V=0; V<_Vfpos.size(); V++){
                for(int f=0; f<_Vfpos[V].size(); f++) {
                    Matrix<Point3>& cmpos = _Vfpos[V][f];
                    for(int j=0; j<HLF; j++) {
                        for(int i=0; i<HLF; i++) {
                            if(mode == GL_SELECT){
                                glLoadName(count);
                            }
                            glBegin(GL_QUADS);
                            glVertex3dv(cmpos(i,j+1));
                            glVertex3dv(cmpos(i,j));
                            glVertex3dv(cmpos(i+1,j));
                            glVertex3dv(cmpos(i+1,j+1));
                            glEnd();
                            count++;
                        }
                    }
                }
            }
            //    }else{
            //        int a = base*base;
            //        int b = base;
            //        for(int ii=0; ii<nfaces; ii++){
            //            Matrix<Point3>& cmpos = _Vftpos[ii];
            //            int mmax = cmpos.m();
            //            int nmax = cmpos.n();
            //            for(int j=0; j<nmax-1; j++) {
            //                for(int i=0; i<mmax-1; i++) {
            //                    if(mode == GL_SELECT){
            //                        count  = a*ii + b*i + j;
            //                        glLoadName(count);
            //                    }
            //                    glBegin(GL_TRIANGLES);
            //
            //                    glVertex3dv(cmpos(i,j+1));
            //                    glVertex3dv(cmpos(i,j));
            //                    glVertex3dv(cmpos(i+1,j));
            //                    //       glVertex3dv(cmpos(i+1,j+1));
            //
            //                    glEnd();
            //                }
            //            }
            //        }
    }
    }


    /*  processHits() prints out the contents of the 
     *  selection array.
     */
    GLuint FaceMapSurfObj::processHits(GLint hits, GLuint buffer[])
    {
        //GLuint names, *ptr;
        GLuint *ptr;
        GLuint ret = 0;

        //  printf("hits = %d\n", hits);
        ptr = (GLuint *) buffer;

        double zmin = 1e6;
        for (int i = 0; i < hits; i++) {  /* for each hit  */
            //names = *ptr;
            ptr++;
            double z1 = *ptr/0x7fffffff; ptr++;
            double z2 = *ptr/0x7fffffff; ptr++;
            if(z2<zmin || z1<zmin){
                zmin = z1;
                ret = *ptr;
            }
            ptr++;
        }
        //  cout << ret << endl;
        return ret;
    }

    void FaceMapSurfObj::setPerspectiveParam(float f, float a, float zn, float zf){
        fovy = f;
        aspect = a;
        znear = zn;
        zfar = zf;
    }
