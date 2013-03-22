/* Some useful links
 * http://www.opencascade.org/org/forum/thread_20775/?forum=3
 * file:///home/fabian/util/opencascade/OpenCASCADE654/doc/ReferenceDocumentation/html/class_topo_d_s___vertex.html#details
 * http://www.opencascade.org/org/forum/thread_16213/?forum=3
 */

#include <iostream>

#include "IGESControl_Controller.hxx"
#include "BRepTools.hxx"
#include "BRep_Tool.hxx"
#include "Geom_Surface.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS_Wire.hxx"
#include "TopoDS_Vertex.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS.hxx"
#include "gp_Pnt.hxx"

#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepAlgoAPI_Common.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepTools_WireExplorer.hxx"

#include "IGESControl_Writer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"

#include "occtFunc.h"

using namespace std;

void commonFace(double XYZL[], double XYZR[], double XYZCommon[], double &AR, double XYZF[])
{
    //points of first face
    gp_Pnt P1F1(XYZL[0],XYZL[1],XYZL[2]);
    gp_Pnt P2F1(XYZL[3],XYZL[4],XYZL[5]);
    gp_Pnt P3F1(XYZL[6],XYZL[7],XYZL[8]);
    gp_Pnt P4F1(XYZL[9],XYZL[10],XYZL[11]);
    TopoDS_Edge aEdge1F1 = BRepBuilderAPI_MakeEdge(P1F1 , P2F1);
    TopoDS_Edge aEdge2F1 = BRepBuilderAPI_MakeEdge(P2F1 , P4F1);
    TopoDS_Edge aEdge3F1 = BRepBuilderAPI_MakeEdge(P1F1 , P3F1);
    TopoDS_Edge aEdge4F1 = BRepBuilderAPI_MakeEdge(P3F1 , P4F1);
    TopoDS_Wire aWireF1 = BRepBuilderAPI_MakeWire(aEdge1F1 , aEdge2F1 , aEdge3F1, aEdge4F1);
    TopoDS_Face aFaceF1 = BRepBuilderAPI_MakeFace(aWireF1);
    TopoDS_Shape S1 = aFaceF1;

    //points of second face
    gp_Pnt P1F2(XYZR[0],XYZR[1],XYZR[2]);
    gp_Pnt P2F2(XYZR[3],XYZR[4],XYZR[5]);
    gp_Pnt P3F2(XYZR[6],XYZR[7],XYZR[8]);
    gp_Pnt P4F2(XYZR[9],XYZR[10],XYZR[11]);
    TopoDS_Edge aEdge1F2 = BRepBuilderAPI_MakeEdge(P1F2 , P2F2);
    TopoDS_Edge aEdge2F2 = BRepBuilderAPI_MakeEdge(P2F2 , P4F2);
    TopoDS_Edge aEdge3F2 = BRepBuilderAPI_MakeEdge(P1F2 , P3F2);
    TopoDS_Edge aEdge4F2 = BRepBuilderAPI_MakeEdge(P3F2 , P4F2);
    TopoDS_Wire aWireF2 = BRepBuilderAPI_MakeWire(aEdge1F2 , aEdge2F2 , aEdge3F2, aEdge4F2);
    TopoDS_Face aFaceF2 = BRepBuilderAPI_MakeFace(aWireF2);
    TopoDS_Shape S2 = aFaceF2;
    
    //calculation of common part
    TopoDS_Shape SC = BRepAlgoAPI_Common(S1,S2);
    

    GProp_GProps sprops;
    BRepGProp::SurfaceProperties(SC, sprops);
    AR=sprops.Mass();
    gp_Pnt comPt = sprops.CentreOfMass();
    TopoDS_Vertex ver_comPt = BRepBuilderAPI_MakeVertex(comPt);
    comPt.Coord(XYZF[0],XYZF[1],XYZF[2]);

    int ii = 0;
    int i_min = 0;
    double x_min, y_min, z_min;
    double XYZC[12];

    if (AR > 1.e-10)
    {
	TopExp_Explorer expShape;
	for (expShape.Init(SC, TopAbs_WIRE);expShape.More();expShape.Next())
	{
	    TopoDS_Wire W = TopoDS::Wire(expShape.Current());
	    BRepTools_WireExplorer ExpWire;
	    for (ExpWire.Init(W); ExpWire.More(); ExpWire.Next())
            {
                gp_Pnt P = BRep_Tool::Pnt(ExpWire.CurrentVertex());
                if ((P.X() < x_min && P.Y() < y_min && P.Z() < z_min) || ii == 0)
                {
                    i_min = ii;
                    x_min = P.X();
                    y_min = P.Y();
                    z_min = P.Z();
                }
                XYZC[ii*3]=P.X();
                XYZC[ii*3+1]=P.Y();
                XYZC[ii*3+2]=P.Z();
                ii++;
            }
        }

        int i_pointer=i_min;
        for (int i=0; i<4; i++)
        {
            if (i_pointer > 3) i_pointer = 0;
            XYZCommon[i*3]=XYZC[i_pointer*3];
            XYZCommon[i*3+1]=XYZC[i_pointer*3+1];
            XYZCommon[i*3+2]=XYZC[i_pointer*3+2];
            i_pointer++;
        }

    }
    else AR=-1.0;
    return;
}
