#ifndef INCLUDE_LOOP_HPP
#define INCLUDE_LOOP_HPP

class loop
{
public:

  // constructor and destructor
  loop(int a_region, int h_region);
  ~loop();

  // setters
  void set_a_bound(bool a_bound);
  void set_h_bound(bool h_bound);
  void set_a(int a);
  void set_h(int h);
  void set_d(int d);
  void set_a_region(int a_region);
  void set_h_region(int h_region);

  // getters
  bool get_a_bound();
  bool get_h_bound();
  int get_a();
  int get_h();
  int get_d();
  int get_a_region();
  int get_h_region();
  
private:

  bool a_bound, h_bound; // binding state of anchor and hinge
  int a, h; // indices of anchor and hinge
  int d; // direction of loop extrusion
  int a_region; // index of region containing anchor
  int h_region; // index of region containing hinge


};

#endif
