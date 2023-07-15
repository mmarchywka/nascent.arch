#ifndef MJM_LIBMESH_ELEM_H__
#define MJM_LIBMESH_ELEM_H__
#include "mjm_globals.h"
#include "mjm_geometry.h"
#include "mjm_libmesh.h"
//#include "mjm_shapes.h"
//#include "mjm_cursing_list.h"
#include "mjm_rational.h"

#include <vector>
//#include <list>
#include <sstream>
#include <string>
#include "libmesh/face_quad4.h"






class mjm_elem : public libMesh::Quad4
{
typedef mjm_elem Myt;
typedef mjm_rational RatTy;
public:
mjm_elem():bi(0) 
{


}

// with one bc node, just apply if new node has SAME x and y ???
// with two nodes, same x or y depending on what is shared.
int existings(BoundaryInfo & bi )
{
	xx.clear(); yy.clear(); xy.clear(); yx.clear();
	int nbct=0;
	std::vector<const Node *> bvns;
	const Elem * that=this->top_parent();
  	const IdxTy nnodes=(*that).n_nodes();
    for (IdxTy i=0; i<nnodes; ++i)
    {
        const Node * nno=(*that).get_node(i);
        const IdxTy nbc=bi.n_boundary_ids(nno);
        if (nbc!=0) { bvns.push_back(nno); ++nbct; }  // m[nno]=i;
		MM_ERR(" for nodes at "<<(*nno)(0)<<" "<<(*nno)(1)<<" "<<nbc<<" nbct so far "<<nbct)
    } // for i 
	if (nbct==0) return 0;
	// find criteria for new bc in children 
	if (nbct==1) 
		//{ xy.push_back((*bvns[0])(0)); yx.push_back((*bvns[0])(1));return nbct; }
		{ xy.push_back(RatTy().set((*bvns[0])(0))); yx.push_back(RatTy().set((*bvns[0])(1)));return nbct; }
		//{ xx.push_back(RatTy().set((*bvns[0])(0))); yy.push_back(RatTy().set((*bvns[0])(1)));return nbct; }
	// they should have shared component, just ignore if not 
	if (nbct==2)
	{
		//if ((*bvns[0])(0)==(*bvns[1])(0)) xx.push_back((*bvns[0])(0));
		if (RatTy().set((*bvns[0])(0)).equals((*bvns[1])(0))) 
			xx.push_back(RatTy().set((*bvns[0])(0)));
		//if ((*bvns[0])(1)==(*bvns[1])(1)) yy.push_back((*bvns[0])(1));
		if (RatTy().set((*bvns[0])(1)).equals(RatTy().set((*bvns[1])(1)))) 
			yy.push_back(RatTy().set((*bvns[0])(1)));
	}	
	// inside corder, need either of x or y to be a hit
	if (nbct==3)
	{
	/* 	if ((*bvns[0])(0)==(*bvns[1])(0)) xx.push_back((*bvns[0])(0));
		if ((*bvns[0])(1)==(*bvns[1])(1)) yy.push_back((*bvns[0])(1));
		if ((*bvns[0])(0)==(*bvns[2])(0)) xx.push_back((*bvns[0])(0));
		if ((*bvns[0])(1)==(*bvns[2])(1)) yy.push_back((*bvns[0])(1));
		if ((*bvns[2])(0)==(*bvns[1])(0)) xx.push_back((*bvns[1])(0));
		if ((*bvns[2])(1)==(*bvns[1])(1)) yy.push_back((*bvns[1])(1));
*/

		setbvns(bvns,xx,0,1,0);
		setbvns(bvns,yy,0,1,1);

		setbvns(bvns,xx,0,2,0);
		setbvns(bvns,yy,0,2,1);
		
		setbvns(bvns,xx,2,1,0);
		setbvns(bvns,yy,2,1,1);


	}	

	return nbct;
}

template <class Tx, class Ty> void setbvns(Tx & bvns, Ty & xxx, const IdxTy i, const IdxTy j, const IdxTy k)
{
	//if ((*bvns[0])(0)==(*bvns[1])(0)) xx.push_back((*bvns[0])(0));
	if (RatTy().set((*bvns[i])(k)).equals(RatTy().set((*bvns[j])(k))))
		 xxx.push_back(RatTy().set((*bvns[i])(k)));
}

bool make_bc(const Node * newnode)
{
const double x=(*newnode)(0);
const double y=(*newnode)(1);

RatTy xr=RatTy().set(x);
RatTy yr=RatTy().set(y);

MM_ERR(" look for  find bc for "<<xr.string()<<" "<<yr.string()<<" in " )

for (IdxTy i=0; i<xy.size(); ++i)
{
	//if ((x==xy[i])&&(y==yx[i])) return true; 
	if ((xr.equals(xy[i]))&&(yr.equals(yx[i]))) return true; 
} //i 
// MM_ERR("child start ")
for (IdxTy i=0; i<xx.size(); ++i) { // { MM_ERR("xxxxx "<<x<<" "<<xx[i]) } 
//if (x==xx[i]) return true; }
if (xr.equals(xx[i])) return true; }
for (IdxTy i=0; i<yy.size(); ++i) { // { MM_ERR("yyyyy "<<y<<" "<<yy[i]) }
  //if (y==yy[i]) return true; }
  if (yr.equals(yy[i])) return true; }
// MM_ERR("child ends false ")

MM_ERR(" fail to find bc for "<<xr.string()<<" "<<yr.string()<<" in " )
for (IdxTy i=0; i<xy.size(); ++i) MM_ERR(" xy i "<<i<<" "<<xy[i].string()<<" "<<yx[i].string())
for (IdxTy i=0; i<xx.size(); ++i) MM_ERR(" xx i "<<i<<" "<<xx[i].string())
for (IdxTy i=0; i<yy.size(); ++i) MM_ERR(" yy i "<<i<<" "<<yy[i].string())


return false; 
}

 virtual void refine (libMesh::MeshRefinement & mesh_refinement)
{ 
// Elem::refine(mesh_refinement); 
MM_ERR(" refine called")
  libmesh_assert_equal_to (this->refinement_flag(), Elem::REFINE);
  libmesh_assert (this->active());
	int do_bc=0;
	if (bi!=0) do_bc=existings(*bi);
  // Create my children if necessary
  if (!_children)
    {
      _children = new Elem *[this->n_children()];

      unsigned int parent_p_level = this->p_level();
      for (unsigned int c=0; c<this->n_children(); c++)
        {
          _children[c] = Elem::build(this->type(), this).release();
          //Elem * current_child = this->child(c);
          Myt * current_child = (Myt*)this->child(c);
			current_child->bi=bi;
          current_child->set_refinement_flag(Elem::JUST_REFINED);
          current_child->set_p_level(parent_p_level);
          current_child->set_p_refinement_flag(this->p_refinement_flag());

          for (unsigned int nc=0; nc<current_child->n_nodes(); nc++)
            {
              Node * node =
                mesh_refinement.add_node(*this, c, nc,
                                         current_child->processor_id());
              node->set_n_systems (this->n_systems());
              current_child->set_node(nc) = node;
			  if (do_bc!=0){
				 //if  (make_bc(node)||(((*node)(1)==.7)))
				 if  (make_bc(node))
					{
				// for now jus add a dummy value 
MM_ERR(" addig bc at "<<(*node)(0)<<" "<<(*node)(1) <<" "<<((Myt*)top_parent())->kids.size())
				(*bi).add_node(node,0); 
				((Myt*)top_parent())->kids.push_back(node);
// this does not fudding allow children to be added what the fudd??????
//				for (IdxTy i=0; i<current_child->n_sides(); ++i)
//				{ bi->add_side(current_child,i, bi->boundary_id(this->top_parent(),i)); }	

					}}


            }

          mesh_refinement.add_elem (current_child);
          current_child->set_n_systems(this->n_systems());
        }
    }
  else

  {
      unsigned int parent_p_level = this->p_level();
      for (unsigned int c=0; c<this->n_children(); c++)
        {
          Elem * current_child = this->child(c);
// wtf would these already exist???

          for (unsigned int nc=0; nc<current_child->n_nodes(); nc++)
			{
              Node * node = current_child->get_node(nc);
			  if (do_bc!=0){
				 if  (make_bc(node))
					{
				// for now jus add a dummy value 
MM_ERR(" addig  exist node bc at "<<(*node)(0)<<" "<<(*node)(1) <<" "<<((Myt*)top_parent())->kids.size())
				(*bi).add_node(node,0); 
				((Myt*)top_parent())->kids.push_back(node);
				}}

			} // nc


          libmesh_assert(current_child->subactive());
          current_child->set_refinement_flag(Elem::JUST_REFINED);
          current_child->set_p_level(parent_p_level);
          current_child->set_p_refinement_flag(this->p_refinement_flag());
        }
    }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::INACTIVE);

  // Leave the p refinement flag set - we will need that later to get
  // projection operations correct
  // this->set_p_refinement_flag(Elem::INACTIVE);

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      libmesh_assert_equal_to (this->child(c)->parent(), this);
      libmesh_assert(this->child(c)->active());
    }
  libmesh_assert (this->ancestor());


} // method  

//private:
BoundaryInfo * bi;
//typedef std::vector<double> Hit;
typedef std::vector<RatTy> Hit;
Hit xx,yy,xy,yx;
std::vector<const Node * > kids;
}; 



#endif

