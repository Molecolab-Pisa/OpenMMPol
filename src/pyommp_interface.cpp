#include "openmmpol.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <string>

namespace py = pybind11;
typedef typename py::array_t<int, py::array::c_style | py::array::forcecast> py_ciarray;
typedef typename py::array_t<double, py::array::c_style | py::array::forcecast> py_cdarray;

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

        void set_external_field(py_cdarray ext_field, 
                                bool nomm = false, 
                                int32_t solver = OMMP_SOLVER_DEFAULT){
            if(ext_field.ndim() != 2 || 
               ext_field.shape(0) != get_pol_atoms() ||
               ext_field.shape(1) != 3){
                throw py::value_error("ext_field should be shaped [pol_atoms, 3]");
            }

            if(! nomm)
                ommp_set_external_field(handler, ext_field.data(), solver);
            else
                ommp_set_external_field_nomm(handler, ext_field.data(), solver);
            return ;
        }



    private: 
        void *handler;
};

PYBIND11_MODULE(pyopenmmpol, m){
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
        .def("set_external_field", 
             &OMMPSystem::set_external_field,
             "Set the external electric field and solves the linear system.",
             py::arg("external_field"),
             py::arg("nomm") = false,
             py::arg("solver") = OMMP_SOLVER_DEFAULT)
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
