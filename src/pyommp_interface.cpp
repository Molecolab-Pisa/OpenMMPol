#include "openmmpol.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>

namespace py = pybind11;

std::map<std::string, int32_t> solvers{
    {"default", OMMP_SOLVER_DEFAULT},
    {"conjugate gradient", OMMP_SOLVER_CG},
    {"cg", OMMP_SOLVER_CG},
    {"inversion", OMMP_SOLVER_INVERSION},
    {"diis", OMMP_SOLVER_DIIS}
};

typedef typename py::array_t<int, py::array::c_style | py::array::forcecast> py_ciarray;
typedef typename py::array_t<double, py::array::c_style | py::array::forcecast> py_cdarray;

std::map<std::string, py_cdarray> charges_elec_prop(py_cdarray csource,
                                                    py_cdarray qsource,
                                                    py_cdarray ctarget,
                                                    bool V, bool E, bool Egrd,
                                                    bool HE){
    // Interface to ommputils_q_elec_prop
    if(ctarget.ndim() != 2 || 
       ctarget.shape(1) != 3){
        throw py::value_error("ctarget should be shaped [:, 3]");
    }
    int nt = ctarget.shape(0);

    if(csource.ndim() != 2 || 
       csource.shape(1) != 3){
        throw py::value_error("csource should be shaped [:, 3]");
    }
    int ns = csource.shape(0);
    
    if(qsource.ndim() != 1 || 
       qsource.shape(0) != ns){
        throw py::value_error("qsource should be shaped [csource.shape[0]]");
    }

    double *Vmem, *Emem, *Egrdmem, *HEmem;
    py_cdarray Varr, Earr, Egrdarr, HEarr;


    if(V)
        Vmem = new double[nt];
    else
        Vmem = NULL;
    if(E)
        Emem = new double[3*nt];
    else
        Emem = NULL;
    if(Egrd)
        Egrdmem = new double[6*nt];
    else
        Egrdmem = NULL;
    if(HE)
        HEmem = new double[10*nt];
    else
        HEmem = NULL;
    
    ommputils_q_elec_prop(ns, nt, csource.data(), qsource.data(), 
                          ctarget.data(), Vmem, Emem, Egrdmem, HEmem);
  
    if(V){
        py::buffer_info buf_V(Vmem, sizeof(double),
                          py::format_descriptor<double>::format(),
                          1,
                          {nt},
                          {sizeof(double)});
        Varr = py_cdarray(buf_V);
    }
    else{
        Varr = py::none();
    }
    if(E){
        py::buffer_info buf_E(Emem, sizeof(double),
                          py::format_descriptor<double>::format(),
                          2,
                          {nt,3},
                          {3*sizeof(double), sizeof(double)});
        Earr = py_cdarray(buf_E);
    }
    else{
        Earr = py::none();
    }
    if(Egrd){
        py::buffer_info buf_Egrd(Egrdmem, sizeof(double),
                             py::format_descriptor<double>::format(),
                             2,
                             {nt, 6},
                             {6*sizeof(double), sizeof(double)});
        Egrdarr = py_cdarray(buf_Egrd);
    }
    else{
        Egrdarr = py::none();
    }
    if(HE){
        py::buffer_info buf_HE(HEmem, sizeof(double),
                          py::format_descriptor<double>::format(),
                          2,
                          {nt, 10},
                          {10*sizeof(double), sizeof(double)});
        HEarr = py_cdarray(buf_HE);
    }
    else{
        HEarr = py::none();
    }

    std::map<std::string, py_cdarray> res{
        {"V", Varr},
        {"E", Earr},
        {"Egrad", Egrdarr},
        {"EHess", HEarr},
    };
    return res;
}

class OMMPSystem{
    public:
        OMMPSystem(){
            handler = nullptr;
        }

        OMMPSystem(std::string mmp_filename){
            handler = ommp_init_mmp(mmp_filename.c_str());
        }
        
        OMMPSystem(std::string xyz_filename, std::string prm_filename){
            handler = ommp_init_xyz(xyz_filename.c_str(), prm_filename.c_str());
        }

        ~OMMPSystem(){
            if(is_init()) ommp_terminate(handler);
            handler = nullptr;
            return;
        }
        
        bool is_init(){
            return handler != nullptr;
        }

        void print_summary(std::string outfile = ""){
            if(outfile.empty())
                ommp_print_summary(handler);
            else
                ommp_print_summary_to_file(handler, 
                                           outfile.c_str());
            return;
        }

        void save_mmp_file(std::string outfile, int version = 3){
            ommp_save_mmp(handler, outfile.c_str(), version);
            return;
        }

        int get_n_ipd(){
            return ommp_get_n_ipd(handler);
        }

        int get_pol_atoms(){
            return ommp_get_pol_atoms(handler);
        }

        int get_mm_atoms(){
            return ommp_get_mm_atoms(handler);
        }

        int get_ld_cart(){
            return ommp_get_ld_cart(handler);
        }

        bool is_amoeba(){
            return ommp_ff_is_amoeba(handler);
        }

        py_cdarray get_ipd(){
            double *mem = ommp_get_ipd(handler);
            py::buffer_info bufinfo(mem, sizeof(double),
                                 py::format_descriptor<double>::format(),
                                 3,
                                 {get_n_ipd(), get_pol_atoms(), 3},
                                 {get_pol_atoms()*3*sizeof(double), 3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray get_cmm(){
            double *mem = ommp_get_cmm(handler);
            py::buffer_info bufinfo(mem, sizeof(double),
                                 py::format_descriptor<double>::format(),
                                 2,
                                 {get_mm_atoms(), 3},
                                 {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }

        py_cdarray get_cpol(){
            double *mem = ommp_get_cpol(handler);
            py::buffer_info bufinfo(mem, sizeof(double),
                                 py::format_descriptor<double>::format(),
                                 2,
                                 {get_pol_atoms(), 3},
                                 {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray get_q(){
            double *mem = ommp_get_q(handler);
            py::buffer_info bufinfo(mem, sizeof(double),
                                 py::format_descriptor<double>::format(),
                                 2,
                                 {get_mm_atoms(), get_ld_cart()},
                                 {get_ld_cart()*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray get_static_charges(){
            double *mem = ommp_get_q(handler);
            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    1,
                                    {get_mm_atoms()},
                                    {get_ld_cart()*sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray get_static_dipoles(){
            if(! is_amoeba()) return py::none();

            double *mem = ommp_get_q(handler);
            py::buffer_info bufinfo(&(mem[1]), sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {get_mm_atoms(), 3},
                                    {get_ld_cart()*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }

        py_cdarray get_static_quadrupoles(){
            if(! is_amoeba()) return py::none();
            
            double *mem = ommp_get_q(handler);
            py::buffer_info bufinfo(&(mem[4]), sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {get_mm_atoms(), 6},
                                    {get_ld_cart()*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_ciarray get_polar_mm(){
            int32_t *mem = ommp_get_polar_mm(handler);
            py::buffer_info bufinfo(mem, sizeof(int),
                                    py::format_descriptor<int>::format(),
                                    1,
                                    {get_pol_atoms()},
                                    {sizeof(int)});
            return py_cdarray(bufinfo);
        }

        py_cdarray mm_potential_at_external(py_cdarray cext){
            if(cext.ndim() != 2 || 
               cext.shape(1) != 3){
                throw py::value_error("cext should be shaped [:, 3]");
            }
            double *mem = new double[cext.shape(0)];
            for(int i=0; i < cext.shape(0); i++)
                mem[i] = 0.0;

            ommp_potential_mm2ext(handler, cext.shape(0), cext.data(), mem);
            
            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    1,
                                    {cext.shape(0)},
                                    {sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray mmpol_field_at_external(py_cdarray cext){
            if(cext.ndim() != 2 || 
               cext.shape(1) != 3){
                throw py::value_error("cext should be shaped [:, 3]");
            }
            
            int n = cext.shape(0);
            double *mem = new double[n*3];
            for(int i=0; i < n*3; i++)
                mem[i] = 0.0;

            ommp_field_mmpol2ext(handler, cext.shape(0), cext.data(), mem);

            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {n, 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray mm_field_at_external(py_cdarray cext){
            if(cext.ndim() != 2 || 
               cext.shape(1) != 3){
                throw py::value_error("cext should be shaped [:, 3]");
            }
            
            int n = cext.shape(0);
            double *mem = new double[n*3];
            for(int i=0; i < n*3; i++)
                mem[i] = 0.0;

            ommp_field_mm2ext(handler, cext.shape(0), cext.data(), mem);

            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {n, 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray pol_field_at_external(py_cdarray cext){
            if(cext.ndim() != 2 || 
               cext.shape(1) != 3){
                throw py::value_error("cext should be shaped [:, 3]");
            }
            
            int n = cext.shape(0);
            double *mem = new double[n*3];
            for(int i=0; i < n*3; i++)
                mem[i] = 0.0;

            ommp_field_pol2ext(handler, cext.shape(0), cext.data(), mem);

            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {n, 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }

        void set_external_field(py_cdarray ext_field, 
                                bool nomm = false, 
                                std::string solver = "default"){
            if(ext_field.ndim() != 2 || 
               ext_field.shape(0) != get_pol_atoms() ||
               ext_field.shape(1) != 3){
                throw py::value_error("ext_field should be shaped [pol_atoms, 3]");
            }

            if(solvers.find(solver) == solvers.end()){
                throw py::value_error("Selected solver is not available!");
            }

            if(! nomm)
                ommp_set_external_field(handler, ext_field.data(), solvers[solver]);
            else
                ommp_set_external_field_nomm(handler, ext_field.data(), solvers[solver]);
            return ;
        }
        
        void update_coordinates(py_cdarray c){ 
            if(c.ndim() != 2 || 
               c.shape(0) != get_mm_atoms() ||
               c.shape(1) != 3){
                throw py::value_error("ext_field should be shaped [mm_atoms, 3]");
            }

            ommp_update_coordinates(handler, c.data());
            return ;
        }

        void numerical_grad(double (OMMPSystem::*ene_f)(void), double *g, double dd=1e-5){
            double *new_c = new double[get_mm_atoms()*3];
            double *_cmm = ommp_get_cmm(handler);
            double tmp;

            for(int i=0; i < get_mm_atoms()*3; i++)
                new_c[i] = _cmm[i];

            for(int i=0;  i < get_mm_atoms()*3; i++){
                new_c[i] += dd;
                ommp_update_coordinates(handler, new_c);
                tmp = (this->*ene_f)();

                new_c[i] -= 2*dd;
                ommp_update_coordinates(handler, new_c);
                g[i] += (tmp-(this->*ene_f)()) / (2*dd);

                new_c[i] += dd;
                ommp_update_coordinates(handler, new_c);
            }
        }

        py_cdarray get_rotation_grad(py_cdarray E, py_cdarray Egrd){
            double *mem = new double[get_mm_atoms()*3];
            for(int i=0; i < get_mm_atoms()*3; i++)
                mem[i] = 0.0;

            ommp_rotation_geomgrad(handler, E.data(), Egrd.data(), mem);
            
            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {get_mm_atoms(), 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray get_fixedelec_grad(bool numerical=false){
            double *mem = new double[get_mm_atoms()*3];
            for(int i=0; i < get_mm_atoms()*3; i++)
                mem[i] = 0.0;

            if(numerical)
                numerical_grad(&OMMPSystem::get_fixedelec_energy, mem);
            else
                ommp_fixedelec_geomgrad(handler, mem);
            
            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {get_mm_atoms(), 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }

        py_cdarray get_polelec_grad(bool numerical=false){
            double *mem = new double[get_mm_atoms()*3];
            for(int i=0; i < get_mm_atoms()*3; i++)
                mem[i] = 0.0;

            if(numerical)
                numerical_grad(&OMMPSystem::get_polelec_energy, mem);
            else
                ommp_polelec_geomgrad(handler, mem);
            
            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {get_mm_atoms(), 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        py_cdarray get_vdw_grad(bool numerical=false){
            double *mem = new double[get_mm_atoms()*3];
            for(int i=0; i < get_mm_atoms()*3; i++)
                mem[i] = 0.0;

            if(numerical)
                numerical_grad(&OMMPSystem::get_vdw_energy, mem);
            else
                ommp_vdw_geomgrad(handler, mem);
            
            py::buffer_info bufinfo(mem, sizeof(double),
                                    py::format_descriptor<double>::format(),
                                    2,
                                    {get_mm_atoms(), 3},
                                    {3*sizeof(double), sizeof(double)});
            return py_cdarray(bufinfo);
        }
        
        double get_bond_energy(void){
            return ommp_get_bond_energy(handler);
        }

        double get_angle_energy(void){
            return ommp_get_angle_energy(handler);
        }

        double get_angtor_energy(void){
            return ommp_get_angtor_energy(handler);
        }

        double get_strtor_energy(void){
            return ommp_get_strtor_energy(handler);
        }

        double get_strbnd_energy(void){
            return ommp_get_strbnd_energy(handler);
        }

        double get_opb_energy(void){
            return ommp_get_opb_energy(handler);
        }

        double get_pitors_energy(void){
            return ommp_get_pitors_energy(handler);
        }

        double get_torsion_energy(void){
            return ommp_get_torsion_energy(handler);
        }

        double get_tortor_energy(void){
            return ommp_get_tortor_energy(handler);
        }

        double get_urey_energy(void){
            return ommp_get_urey_energy(handler);
        }

        double get_polelec_energy(void){
            return ommp_get_polelec_energy(handler);
        }

        double get_fixedelec_energy(void){
            return ommp_get_fixedelec_energy(handler);
        }

        double get_vdw_energy(void){
            return ommp_get_vdw_energy(handler);
        }

    private: 
        void *handler;
};

void set_verbose(int32_t v){
    ommp_set_verbose(v);
}

PYBIND11_MODULE(pyopenmmpol, m){
    m.def("set_verbose", &set_verbose);
    m.attr("available_solvers") = solvers;
    m.def("charges_elec_prop", &charges_elec_prop, 
          "Compute electrostatic properties of a set of charges at arbitrary coordinates.");

    py::class_<OMMPSystem, std::shared_ptr<OMMPSystem>>(m, "OMMPSystem", "System of OMMP library.")
        .def(py::init<std::string>(), 
             "pyOpenMMPol creator, takes the path to a .mmp file as input.", 
             py::arg("mmp_filename"))
        .def(py::init<std::string, std::string>(), 
             "pyOpenMMPol creator, takes the path to a Tinker .xyz and .prm files as input.", 
             py::arg("xyz_filename"), 
             py::arg("prm_filename"))
        .def("print_summary", 
             &OMMPSystem::print_summary, 
             "Output a summary of loaded quantites, if outfile is specified, it is printed on file.", 
             py::arg("outfile") = "")
        .def("save_mmp_file", 
             &OMMPSystem::save_mmp_file, 
             "Save the data of OMMPSystem into a .mmp file (note that some informations are dropped in the process).", 
             py::arg("outfile"), 
             py::arg("version") = 3)
        .def("mm_potential_at_external",
             &OMMPSystem::mm_potential_at_external,
             "Compute the electrostatic potential of the static part of the system to arbitrary coordinates.",
             py::arg("coord"))
        .def("pol_field_at_external",
             &OMMPSystem::pol_field_at_external,
             "Compute the electrostatic field of the polarizable part of the system to arbitrary coordinates.",
             py::arg("cext"))
        .def("mm_field_at_external",
             &OMMPSystem::mm_field_at_external,
             "Compute the electrostatic field of the static part of the system to arbitrary coordinates.",
             py::arg("cext"))
        .def("mmpol_field_at_external",
             &OMMPSystem::mmpol_field_at_external,
             "Compute the electrostatic field of the static and polarizable part of the system to arbitrary coordinates.",
             py::arg("cext"))

        .def("set_external_field", 
             &OMMPSystem::set_external_field,
             "Set the external electric field and solves the linear system.",
             py::arg("external_field"),
             py::arg("nomm") = false,
             py::arg("solver") = "default")
        
        .def("update_coordinates", 
             &OMMPSystem::update_coordinates,
             "Update system's coordinates.",
             py::arg("coord"))
        .def("do_fixedelec_grad",
             &OMMPSystem::get_fixedelec_grad,
             "Compute the geometrical gradients of fixed electrostatic energy.",
             py::arg("numerical") = false)
        .def("do_polelec_grad",
             &OMMPSystem::get_polelec_grad,
             "Compute the geometrical gradients of polarizable electrostatic energy.",
             py::arg("numerical") = false)
        .def("do_rotation_grad",
             &OMMPSystem::get_rotation_grad,
             "Compute the geometrical gradients on MM atoms due to the rotation of multipoles in an external electric field.",
             py::arg("electric_field"), py::arg("electric_field_gradients"))
        .def("do_vdw_grad",
             &OMMPSystem::get_vdw_grad,
             "Compute the geometrical gradients of Van der Waal energy.",
             py::arg("numerical") = false)
        .def("get_bond_energy", &OMMPSystem::get_bond_energy, "Compute the energy of bond stretching")
        .def("get_angle_energy", &OMMPSystem::get_angle_energy, "Compute the energy of angle bending")
        .def("get_torsion_energy", &OMMPSystem::get_torsion_energy, "Compute the energy of dihedral torsion")
        .def("get_opb_energy", &OMMPSystem::get_opb_energy, "Compute the energy of out-of-plane bending")
        .def("get_urey_energy", &OMMPSystem::get_urey_energy, "Compute the energy of Urey-Bradley terms")
        .def("get_pitors_energy", &OMMPSystem::get_pitors_energy, "Compute the energy of Pi-system torsion terms")
        .def("get_strbnd_energy", &OMMPSystem::get_strbnd_energy, "Compute the energy of stretching-bending couplings")
        .def("get_strtor_energy", &OMMPSystem::get_strtor_energy, "Compute the energy of stretching torsion couplings")
        .def("get_angtor_energy", &OMMPSystem::get_angtor_energy, "Compute the energy of bending torsion couplings")
        .def("get_tortor_energy", &OMMPSystem::get_tortor_energy, "Compute the energy of torsion-torsion cmap")
        
        .def("get_vdw_energy", &OMMPSystem::get_vdw_energy, "Compute the energy of Van der Waals terms")
        .def("get_fixedelec_energy", &OMMPSystem::get_fixedelec_energy, "Compute the energy of fixed electrostatics")
        .def("get_polelec_energy", &OMMPSystem::get_polelec_energy, "Compute the energy of polarizable electrostatics")

        .def_property_readonly("pol_atoms", &OMMPSystem::get_pol_atoms, "Number of polarizable atoms")
        .def_property_readonly("mm_atoms", &OMMPSystem::get_mm_atoms, "Number of atoms")
        .def_property_readonly("_ld_cart", &OMMPSystem::get_ld_cart, "Dimension of static site descriptor; 1 for Wang FF, 10 for Amoeba FF.")
        .def_property_readonly("n_ipd", &OMMPSystem::get_n_ipd, "Number of induced dipoles sets; 1 for Wang FF, 2 for Amoeba FF.")
        .def_property_readonly("is_amoeba", &OMMPSystem::is_amoeba, "Flag for Amoeba FF")
        .def_property_readonly("ipd", &OMMPSystem::get_ipd, "Induced dipoles (read only)")
        .def_property_readonly("cmm", &OMMPSystem::get_cmm, "Coordinates of atoms (read only) [mm_atoms, 3]")
        .def_property_readonly("cpol", &OMMPSystem::get_cpol, "Coordinates of polarizable atoms (read only) [pol_atoms, 3]")
        .def_property_readonly("_q", &OMMPSystem::get_q, "Full descriptor of static sites (read only) [mm_atoms, _ld_cart]")
        .def_property_readonly("static_charges", &OMMPSystem::get_static_charges, "Static charges (read only) [mm_atoms]")
        .def_property_readonly("static_dipoles", &OMMPSystem::get_static_dipoles, "Static dipoles (read only) if is_amoeba [mm_atoms, 3], else None")
        .def_property_readonly("static_quadrupoles", &OMMPSystem::get_static_quadrupoles, "Static dipoles (read only) if is_amoeba [mm_atoms, 6], else None")
        .def_property_readonly("polar_mm", &OMMPSystem::get_polar_mm, "Index of polarizable atoms in atom list [pol_atoms]");

}
