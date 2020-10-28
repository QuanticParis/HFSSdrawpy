import atexit
from copy import copy
import os

import tempfile
import types
import numpy
import signal
import pythoncom
import time
from functools import wraps
from sympy.parsing import sympy_parser
from pint import UnitRegistry
from win32com.client import Dispatch, CDispatch

from ..utils import parse_entry, \
                            val, \
                            LENGTH_UNIT, \
                            Vector
                            #extract_value_unit, \
                            #extract_value_dim, \
                            #rem_unit, \

ureg = UnitRegistry()
Q = ureg.Quantity

BASIS_ORDER = {"Zero Order": 0,
               "First Order": 1,
               "Second Order": 2,
               "Mixed Order": -1}

layer_Default = 10

def increment_name(base, existing):
    if base not in existing:
        return base
    n = 1
    make_name = lambda: base + str(n)
    while make_name() in existing:
        n += 1
    return make_name()

_release_fns = []

def _add_release_fn(fn):
    global _release_fns
    _release_fns.append(fn)
    atexit.register(fn)
    signal.signal(signal.SIGTERM, fn)
    signal.signal(signal.SIGABRT, fn)


def release():
    global _release_fns
    for fn in _release_fns:
        fn()
    time.sleep(0.1)
    refcount = pythoncom._GetInterfaceCount()
    if refcount > 0:
        print("Warning! %d COM references still alive")
        print("HFSS will likely refuse to shut down")


class COMWrapper(object):
    def __init__(self):
        _add_release_fn(self.release)

    def release(self):
        for k, v in self.__dict__.items():
            if isinstance(v, CDispatch):
                setattr(self, k, None)


class HfssPropertyObject(COMWrapper):
    prop_holder = None
    prop_tab = None
    prop_server = None


def make_str_prop(name, prop_tab=None, prop_server=None):
    return make_prop(name, prop_tab=prop_tab, prop_server=prop_server)


def make_int_prop(name, prop_tab=None, prop_server=None):
    return make_prop(name, prop_tab=prop_tab, prop_server=prop_server, prop_args=["MustBeInt:=", True])


def make_float_prop(name, prop_tab=None, prop_server=None):
    return make_prop(name, prop_tab=prop_tab, prop_server=prop_server, prop_args=["MustBeInt:=", False])


def make_prop(name, prop_tab=None, prop_server=None, prop_args=None):
    def set_prop(self, value, prop_tab=prop_tab, prop_server=prop_server,
                 prop_args=prop_args):
        prop_tab = self.prop_tab if prop_tab is None else prop_tab
        prop_server = self.prop_server if prop_server is None else prop_server
        if isinstance(prop_tab, types.FunctionType):
            prop_tab = prop_tab(self)
        if isinstance(prop_server, types.FunctionType):
            prop_server = prop_server(self)
        if prop_args is None:
            prop_args = []
        self.prop_holder.ChangeProperty(
            ["NAME:AllTabs",
             ["NAME:"+prop_tab,
              ["NAME:PropServers", prop_server],
              ["NAME:ChangedProps",
               ["NAME:"+name, "Value:=", value] + prop_args]]])

    def get_prop(self, prop_tab=prop_tab, prop_server=prop_server):
        prop_tab = self.prop_tab if prop_tab is None else prop_tab
        prop_server = self.prop_server if prop_server is None else prop_server
        if isinstance(prop_tab, types.FunctionType):
            prop_tab = prop_tab(self)
        if isinstance(prop_server, types.FunctionType):
            prop_server = prop_server(self)
        return self.prop_holder.GetPropertyValue(prop_tab, prop_server, name)

    return property(get_prop, set_prop)


class HfssApp(COMWrapper):
    def __init__(self):
        super(HfssApp, self).__init__()
        self._app = Dispatch('AnsoftHfss.HfssScriptInterface')

    def get_app_desktop(self):
        return HfssDesktop(self, self._app.GetAppDesktop())


class HfssDesktop(COMWrapper):
    def __init__(self, app, desktop):
        """
        :type app: HfssApp
        :type desktop: Dispatch
        """
        super(HfssDesktop, self).__init__()
        self.parent = app
        self._desktop = desktop

    def close_all_windows(self):
        self._desktop.CloseAllWindows()

    def project_count(self):
        return self._desktop.Count()

    def get_active_project(self):
        return HfssProject(self, self._desktop.GetActiveProject())

    def get_projects(self):
        return [HfssProject(self, p) for p in self._desktop.GetProjects()]

    def get_project_names(self):
        return self._desktop.GetProjectList()

    def get_version(self):
        return self._desktop.GetVersion()

    def new_project(self):
        return HfssProject(self, self._desktop.NewProject())

    def open_project(self, path):
        return HfssProject(self, self._desktop.OpenProject(path))

    def set_active_project(self, name):
        self._desktop.SetActiveProject(name)

    def clear_all_messages(self):
        self._desktop.ClearMessages("", "", 3)

    @property
    def project_directory(self):
        return self._desktop.GetProjectDirectory()

    @project_directory.setter
    def project_directory(self, path):
        self._desktop.SetProjectDirectory(path)

    @property
    def library_directory(self):
        return self._desktop.GetLibraryDirectory()

    @library_directory.setter
    def library_directory(self, path):
        self._desktop.SetLibraryDirectory(path)

    @property
    def temp_directory(self):
        return self._desktop.GetTempDirectory()

    @temp_directory.setter
    def temp_directory(self, path):
        self._desktop.SetTempDirectory(path)


class HfssProject(COMWrapper):
    def __init__(self, desktop, project):
        """
        :type desktop: HfssDesktop
        :type project: Dispatch
        """
        super(HfssProject, self).__init__()
        self.parent = desktop
        self._project = project
        self.name = project.GetName()

    def close(self):
        self._project.Close()

    def make_active(self):
        self.parent.set_active_project(self.name)

    def get_designs(self):
        return [HfssDesign(self, d) for d in self._project.GetDesigns()]

    def save(self, path=None):
        if path is None:
            self._project.Save()
        else:
            self._project.SaveAs(path, True)

    def simulate_all(self):
        self._project.SimulateAll()

    def import_dataset(self, path):
        self._project.ImportDataset(path)

    def rename_design(self, design, rename):
        if design in self.get_designs():
            design.rename_design(design.name, rename)
        else:
            raise ValueError('%s design does not exist' % design.name)

    def duplicate_design(self, target, source):
        src_design = self.get_design(source)
        return src_design.duplicate(name=target)

    def get_variable_names(self):
        return [s for s in self._project.GetVariables()]

    def get_variables(self):
        return {s: self.get_variable_value(s) for s in self._project.GetVariables()}

    def get_variable_value(self, name):
        return self._project.GetVariableValue(name)

    def create_variable(self, name, value):
        self._project.ChangeProperty(
            ["NAME:AllTabs",
             ["NAME:ProjectVariableTab",
              ["NAME:PropServers", "ProjectVariables"],
              ["Name:NewProps",
               ["NAME:" + name,
                "PropType:=", "VariableProp",
                "UserDef:=", True,
                "Value:=", value]]]])

    def set_variable(self, name, value):
#        print(self._project.GetVariables())
        if name not in self._project.GetVariables():
            self.create_variable(name, value)
        else:
            self._project.SetVariableValue(name, value)

    def get_path(self):
        return self._project.GetPath()

    def new_design(self, name, type):
        name = increment_name(name, [d.GetName() for d in self._project.GetDesigns()])
        return HfssDesign(self, self._project.InsertDesign("HFSS", name, type, ""))

    def get_design(self, name):
        return HfssDesign(self, self._project.GetDesign(name))

    def get_active_design(self):
        d = self._project.GetActiveDesign()
        if d is None:
            raise EnvironmentError("No Design Active")
        return HfssDesign(self, d)

    def new_dm_design(self, name):
        return self.new_design(name, "DrivenModal")

    def new_em_design(self, name):
        return self.new_design(name, "Eigenmode")


class HfssDesign(COMWrapper):
    def __init__(self, project, design):
        super(HfssDesign, self).__init__()
        self.parent = project
        self._design = design
        self.name = design.GetName()
        self.solution_type = design.GetSolutionType()
        if design is None:
            return
        self._setup_module = design.GetModule("AnalysisSetup")
        self._solutions = design.GetModule("Solutions")
        self._fields_calc = design.GetModule("FieldsReporter")
        self._output = design.GetModule("OutputVariable")
        self._boundaries = design.GetModule("BoundarySetup")
        self._mesh = design.GetModule("MeshSetup")
        self._reporter = design.GetModule("ReportSetup")
        self._modeler = design.SetActiveEditor("3D Modeler")
        self._optimetrics = design.GetModule("Optimetrics")
        self.modeler = HfssModeler(self, self._modeler, self._boundaries,
                                   self._mesh)
        self.variables = {}

    def rename_design(self, name):
        old_name = self._design.GetName()
        self._design.RenameDesignInstance(old_name, name)

    def copy_to_project(self, project):
        project.make_active()
        project._project.CopyDesign(self.name)
        project._project.Paste()
        return project.get_active_design()

    def duplicate(self, name=None):
        dup = self.copy_to_project(self.parent)
        if name is not None:
            dup.rename_design(name)
        return dup

    def get_setup_names(self):
        return self._setup_module.GetSetups()

    def get_setup(self, name=None):
        """
        :rtype: HfssSetup
        """
        setups = self.get_setup_names()
        if not setups:
            raise EnvironmentError("No Setups Present")
        if name is None:
            name = setups[0]
        elif name not in setups:
            raise EnvironmentError("Setup {} not found: {}".format(name,
                                                                   setups))

        if self.solution_type == "Eigenmode":
            return HfssEMSetup(self, name)
        elif self.solution_type == "DrivenModal":
            return HfssDMSetup(self, name)

    def create_dm_setup(self, freq_ghz=1, name="Setup", max_delta_s=0.1,
                        max_passes=10, min_passes=1, min_converged=1,
                        pct_refinement=30, basis_order=-1):

        name = increment_name(name, self.get_setup_names())
        self._setup_module.InsertSetup(
            "HfssDriven", [
                "NAME:"+name,
                "Frequency:=", str(freq_ghz)+"GHz",
                "MaxDeltaS:=", max_delta_s,
                "MaximumPasses:=", max_passes,
                "MinimumPasses:=", min_passes,
                "MinimumConvergedPasses:=", min_converged,
                "PercentRefinement:=", pct_refinement,
                "IsEnabled:=", True,
                "BasisOrder:=", basis_order
            ])
        return HfssDMSetup(self, name)

    def create_em_setup(self, name="Setup", min_freq_ghz=1, n_modes=1,
                        max_delta_f=0.1, max_passes=10, min_passes=1,
                        min_converged=1, pct_refinement=30, basis_order=-1):

        name = increment_name(name, self.get_setup_names())
        self._setup_module.InsertSetup(
            "HfssEigen", [
                "NAME:"+name,
                "MinimumFrequency:=", str(min_freq_ghz)+"GHz",
                "NumModes:=", n_modes,
                "MaxDeltaFreq:=", max_delta_f,
                "ConvergeOnRealFreq:=", True,
                "MaximumPasses:=", max_passes,
                "MinimumPasses:=", min_passes,
                "MinimumConvergedPasses:=", min_converged,
                "PercentRefinement:=", pct_refinement,
                "IsEnabled:=", True,
                "BasisOrder:=", basis_order
            ])
        return HfssEMSetup(self, name)

    def delete_setup(self, name):
        if name in self.get_setup_names():
            self._setup_module.DeleteSetups(name)

    def get_nominal_variation(self):
        return self._design.GetNominalVariation()

    def create_variable(self, name, value, postprocessing=False):
        if postprocessing==True:
            variableprop = "PostProcessingVariableProp"
        else:
            variableprop = "VariableProp"

        self._design.ChangeProperty(
            ["NAME:AllTabs",
             ["NAME:LocalVariableTab",
              ["NAME:PropServers", "LocalVariables"],
              ["Name:NewProps",
               ["NAME:" + name,
                "PropType:=", variableprop,
                "UserDef:=", True,
                "Value:=", str(value)]]]])

    def set_variable(self, name, value, postprocessing=False):
        if name not in self._design.GetVariables()+self._design.GetPostProcessingVariables():
            self.create_variable(name, value, postprocessing=postprocessing)
        else:
            self._design.SetVariableValue(name, value)

    def get_variable_value(self, name):
        return self._design.GetVariableValue(name)

    def get_variable_names(self):
        return [s for s in self._design.GetVariables()+self._design.GetPostProcessingVariables()]

    def get_variables(self):
        local_variables = self._design.GetVariables()+self._design.GetPostProcessingVariables()
        return {lv : self.get_variable_value(lv) for lv in local_variables}

    def copy_design_variables(self, source_design):
        ''' does not check that variables are all present '''

        # don't care about values
        source_variables = source_design.get_variables()

        for name, value in source_variables.iteritems():
            self.set_variable(name, value)

    def get_excitations(self):
        self._boundaries.GetExcitations()

    def _evaluate_variable_expression(self, expr, units):
        """
        :type expr: str
        :type units: str
        :return: float
        """
        try:
            sexp = sympy_parser.parse_expr(expr)
        except SyntaxError:
            return Q(expr).to(units).magnitude

        sub_exprs = {fs: self.get_variable_value(fs.name) for fs in sexp.free_symbols}
        return float(sexp.subs({fs: self._evaluate_variable_expression(e, units) for fs, e in sub_exprs.items()}))

    def eval_expr(self, expr, units="mm"):
        return str(self._evaluate_variable_expression(expr, units)) + units

    def eval_var_str(self, name, unit=None): # can only parse 2 times for now
        # TODO: parse several times
        try:
            _val = float(eval(str(sympy_parser.parse_expr(str(sympy_parser.parse_expr(str(name), self.variables)), self.variables))))
        except Exception:
            msg = ('Parsed expression contains a string which does '
                             'not correspond to any design variable')
            raise ValueError(msg)
        if unit is not None:
            _val = str(_val)+' '+unit
        return _val

    def Clear_Field_Clac_Stack(self):
        self._fields_calc.CalcStack("Clear")

class HfssSetup(HfssPropertyObject):
    prop_tab = "HfssTab"
    passes = make_int_prop("Passes")
    pct_refinement = make_float_prop("Percent Refinement")
    basis_order = make_str_prop("Basis Order")

    def __init__(self, design, setup):
        """
        :type design: HfssDesign
        :type setup: Dispatch
        """
        super(HfssSetup, self).__init__()
        self.parent = design
        self.prop_holder = design._design
        self._setup_module = design._setup_module
        self._reporter = design._reporter
        self._solutions = design._solutions
        self.name = setup
        self.solution_name = setup + " : LastAdaptive"
        self.prop_server = "AnalysisSetup:" + setup
        self.expression_cache_items = []

    def analyze(self, name=None):
        if name is None:
            name = self.name
        self.parent._design.Analyze(name)

    def insert_sweep(self, start_ghz, stop_ghz, count=None, step_ghz=None,
                     name="Sweep", type="Fast", save_fields=False):
        if (count is None) == (step_ghz is None):
            raise ValueError("Exactly one of 'points' and 'delta' must be specified")
        name = increment_name(name, self.get_sweep_names())
        params = [
            "NAME:"+name,
            "IsEnabled:=", True,
            "StartValue:=", "%fGHz" % start_ghz,
            "StopValue:=", "%fGHz" % stop_ghz,
            "Type:=", type,
            "SaveFields:=", save_fields,
            "ExtrapToDC:=", False,
        ]
        if step_ghz is not None:
            params.extend([
                "SetupType:=", "LinearSetup",
                "StepSize:=", "%fGHz" % step_ghz,
            ])
        else:
            params.extend([
                "SetupType:=", "LinearCount",
                "Count:=", count,
            ])

        self._setup_module.InsertFrequencySweep(self.name, params)
        return HfssFrequencySweep(self, name)

    def delete_sweep(self, name):
        self._setup_module.DeleteSweep(self.name, name)

    def add_fields_convergence_expr(self, expr, pct_delta, phase=0):
        """note: because of hfss idiocy, you must call "commit_convergence_exprs" after adding all exprs"""
        assert isinstance(expr, NamedCalcObject)
        self.expression_cache_items.append(
            ["NAME:CacheItem",
             "Title:=", expr.name+"_conv",
             "Expression:=", expr.name,
             "Intrinsics:=", "Phase='{}deg'".format(phase),
             "IsConvergence:=", True,
             "UseRelativeConvergence:=", 1,
             "MaxConvergenceDelta:=", pct_delta,
             "MaxConvergeValue:=", "0.05",
             "ReportType:=", "Fields",
             ["NAME:ExpressionContext"]])

    def commit_convergence_exprs(self):
        """note: this will eliminate any convergence expressions not added through this interface"""
        args = [
            "NAME:"+self.name,
            ["NAME:ExpressionCache", self.expression_cache_items]
        ]
        self._setup_module.EditSetup(self.name, args)

    def get_sweep_names(self):
        return self._setup_module.GetSweeps(self.name)

    def get_sweep(self, name=None):
        sweeps = self.get_sweep_names()
        if not sweeps:
            raise EnvironmentError("No Sweeps Present")
        if name is None:
            name = sweeps[0]
        elif name not in sweeps:
            raise EnvironmentError("Sweep {} not found in {}".format(name, sweeps))
        return HfssFrequencySweep(self, name)

    def get_convergence(self, variation=""):
        fn = tempfile.mktemp()
        self.parent._design.ExportConvergence(self.name, variation, fn, False)
        return numpy.loadtxt(fn)

    def get_mesh_stats(self, variation=""):
        fn = tempfile.mktemp()
        self.parent._design.ExportMeshStats(self.name, variation, fn, False)
        return numpy.loadtxt(fn)

    def get_profile(self, variation=""):
        fn = tempfile.mktemp()
        self.parent._design.ExportProfile(self.name, variation, fn, False)
        return numpy.loadtxt(fn)

    def get_fields(self):
        return HfssFieldsCalc(self)

class HfssDMSetup(HfssSetup):
    solution_freq = make_float_prop("Solution Freq")
    delta_s = make_float_prop("Delta S")
    solver_type = make_str_prop("Solver Type")

    def setup_link(self, linked_setup):
        '''
            type: linked_setup <HfssSetup>
        '''
        args = ["NAME:" + self.name,
                ["NAME:MeshLink",
                "Project:=", "This Project*",
                "Design:=", linked_setup.parent.name,
                "Soln:=", linked_setup.solution_name,
                self._map_variables_by_name(),
                "ForceSourceToSolve:=", True,
                "PathRelativeTo:=", "TargetProject",
                ],
                ]
        self._setup_module.EditSetup(self.name, args)

    def _map_variables_by_name(self):
        ''' does not check that variables are all present '''
        # don't care about values
        project_variables = self.parent.parent.get_variable_names()
        design_variables = self.parent.get_variable_names()

        # build array
        args = ["NAME:Params",]
        for name in project_variables:
            args.extend([str(name)+":=", str(name)])
        for name in design_variables:
            args.extend([str(name)+":=", str(name)])
        return args

    def get_solutions(self):
        return HfssDMDesignSolutions(self, self.parent._solutions)

class HfssEMSetup(HfssSetup):
    min_freq = make_float_prop("Min Freq")
    n_modes = make_int_prop("Modes")
    delta_f = make_float_prop("Delta F")

    def get_solutions(self):
        return HfssEMDesignSolutions(self, self.parent._solutions)

class HfssDesignSolutions(COMWrapper):
    def __init__(self, setup, solutions):
        '''
        :type setup: HfssSetup
        '''
        super(HfssDesignSolutions, self).__init__()
        self.parent = setup
        self._solutions = solutions

class HfssEMDesignSolutions(HfssDesignSolutions):
    def eigenmodes(self, lv=""):
        fn = tempfile.mktemp()
        self._solutions.ExportEigenmodes(self.parent.solution_name, lv, fn)
        data = numpy.loadtxt(fn, dtype='str')
        if numpy.size(numpy.shape(data)) == 1: # getting around the very annoying fact that
            data = numpy.array([data])         # in Python a 1D array does not have shape (N,1)
        else:                                  # but rather (N,) ....
            pass
        if numpy.size(data[0,:])==6: # checking if values for Q were saved
            kappa_over_2pis = [2*float(ii) for ii in data[:,3]] # eigvalue=(omega-i*kappa/2)/2pi
                                                    # so kappa/2pi = 2*Im(eigvalue)
        else:
            kappa_over_2pis = None

        freqs = [float(ii) for ii in data[:,1]]
        return freqs, kappa_over_2pis

    def set_mode(self, n, phase):
        n_modes = int(self.parent.n_modes)
        self._solutions.EditSources(
            "TotalFields",
            ["NAME:SourceNames", "EigenMode"],
            ["NAME:Modes", n_modes],
            ["NAME:Magnitudes"] + [1 if i + 1 == n else 0 for i in range(n_modes)],
            ["NAME:Phases"] + [phase if i + 1 == n else 0 for i in range(n_modes)],
            ["NAME:Terminated"], ["NAME:Impedances"]
        )

class HfssDMDesignSolutions(HfssDesignSolutions):
    pass

class HfssFrequencySweep(COMWrapper):
    prop_tab = "HfssTab"
    start_freq = make_float_prop("Start")
    stop_freq = make_float_prop("Stop")
    step_size = make_float_prop("Step Size")
    count = make_float_prop("Count")
    sweep_type = make_str_prop("Type")

    def __init__(self, setup, name):
        """
        :type setup: HfssSetup
        :type name: str
        """
        super(HfssFrequencySweep, self).__init__()
        self.parent = setup
        self.name = name
        self.solution_name = self.parent.name + " : " + name
        self.prop_holder = self.parent.prop_holder
        self.prop_server = self.parent.prop_server + ":" + name

    def analyze_sweep(self):
        self.parent.analyze(self.solution_name)

    def get_network_data(self, formats):
        if isinstance(formats, str):
            formats = formats.split(",")
        formats = [f.upper() for f in formats]

        fmts_lists = {'S': [], 'Y': [], 'Z': []}

        for f in formats:
            fmts_lists[f[0]].append((int(f[1]), int(f[2])))

        ret = [None] * len(formats)
        freq = None

        for data_type, list in fmts_lists.items():
            if list:
                fn = tempfile.mktemp()
                self.parent._solutions.ExportNetworkData(
                    [],  self.parent.name + " : " + self.name,
                      2, fn, ["all"], False, 0,
                      data_type, -1, 1, 15
                )
                with open(fn) as f:
                    f.readline()
                    colnames = f.readline().split()
                array = numpy.loadtxt(fn, skiprows=2)
                if freq is None:
                    freq = array[:, 0]
                for i, j in list:
                    real_idx = colnames.index("%s[%d,%d]_Real" % (data_type, i, j))
                    imag_idx = colnames.index("%s[%d,%d]_Imag" % (data_type, i, j))
                    c_arr = array[:, real_idx] + 1j*array[:, imag_idx]
                    ret[formats.index("%s%d%d" % (data_type, i, j))] = c_arr

        return freq, ret

    def create_report(self, name, expr):
        existing = self.parent._reporter.GetAllReportNames()
        name = increment_name(name, existing)
        var_names = self.parent.parent.get_variable_names()
        var_args = sum([["%s:=" % v_name, ["Nominal"]] for v_name in var_names], [])
        self.parent._reporter.CreateReport(
            name, "Modal Solution Data", "Rectangular Plot",
            self.solution_name, ["Domain:=", "Sweep"], ["Freq:=", ["All"]] + var_args,
            ["X Component:=", "Freq", "Y Component:=", [expr]], [])
        return HfssReport(self.parent.parent, name)

    def get_report_arrays(self, expr):
        r = self.create_report("Temp", expr)
        return r.get_arrays()


class HfssReport(COMWrapper):
    def __init__(self, design, name):
        """
        :type design: HfssDesign
        :type name: str
        """
        super(HfssReport, self).__init__()
        self.parent_design = design
        self.name = name

    def export_to_file(self, filename):
        filepath = os.path.abspath(filename)
        self.parent_design._reporter.ExportToFile(self.name, filepath)

    def get_arrays(self):
        fn = tempfile.mktemp(suffix=".csv")
        self.export_to_file(fn)
        return numpy.loadtxt(fn, skiprows=1, delimiter=',').transpose()


class HfssModeler(COMWrapper):
    def __init__(self, design, modeler, boundaries, mesh):
        """
        :type design: HfssDesign
        """
        super(HfssModeler, self).__init__()
        self.parent = design
        self._modeler = modeler
        self._boundaries = boundaries
        self._mesh = mesh

    def _attributes_array(self, name=None, nonmodel=False, color=None,
                          transparency=0.9, material=None, solve_inside = None,
                          coor_sys='Global', **kwargs):
        arr = ["NAME:Attributes", "PartCoordinateSystem:=", coor_sys, "Coordinate System:=", coor_sys]
        if name is not None:
            arr.extend(["Name:=", name])
        if nonmodel:
            arr.extend(["Flags:=", "NonModel"])

        if color is not None:
            arr.extend(["Color:=", "(%d %d %d)" % color])
        if transparency is not None:
            arr.extend(["Transparency:=", transparency])
        if material is not None:
            arr.extend(["MaterialName:=", '\"' + material + '\"'])
        if solve_inside is not None:
            arr.extend(["SolveInside:=", solve_inside])
        return arr

    def _selections_array(self, *names):
        return ["NAME:Selections", "Selections:=", ",".join(names)]

    def assert_name(func):
        @wraps(func)
        def asserted_name(*args, **kwargs):
            name = func(*args, **kwargs)
            msg = 'Failed at generating a name for %s'%name
            assert name == kwargs['name'], msg
            return name
        return asserted_name

    def connect_faces(self, entity1, entity2):
        name = self._modeler.Connect(["NAME:Selections",
                               "Selections:=", ','.join([entity1.name, entity2.name])])
        return name

    def copy(self, entity):
        self._modeler.Copy(["NAME:Selections", "Selections:=", entity.name])
        new_obj = self._modeler.Paste()
        return new_obj[0]

    def create_coor_sys(self, coor_sys='chip', rel_coor=None,
                        ref_name='Global'):
        if not(ref_name in self._modeler.GetCoordinateSystems()):
            raise ValueError('%s reference coordinate system does \
                             not exist'%ref_name)
        else:
            self.set_coor_sys(ref_name)
        if rel_coor is None:
            rel_coor = [[0, 0, 0],  # origin
                        [1, 0, 0],  # new_x
                        [0, 1, 0]]  # new_y
        origin = rel_coor[0]
        new_x = rel_coor[1]
        new_y = rel_coor[2]
        if not(coor_sys in self._modeler.GetCoordinateSystems()):
            # if the coor_sys does not exist : create
            self._modeler.CreateRelativeCS(["NAME:RelativeCSParameters",
            "Mode:="		, "Axis/Position",
            "OriginX:=", str(origin[0]), "OriginY:="	, str(origin[1]), "OriginZ:=", str(origin[2]),
            "XAxisXvec:=", str(new_x[0]), "XAxisYvec:=", str(new_x[1]), "XAxisZvec:=", str(new_x[2]),
            "YAxisXvec:=", str(new_y[0]), "YAxisYvec:="		, str(new_y[1]), "YAxisZvec:=", str(new_y[2])],

            ["NAME:Attributes", "Name:=", coor_sys])
        else:
            # if the coor_sys exists : modify
            self._modeler.ChangeProperty(["NAME:AllTabs",
		["NAME:Geometry3DCSTab",
			["NAME:PropServers", coor_sys],
			["NAME:ChangedProps", ["NAME:Reference CS","Value:=", ref_name],
                                ["NAME:Origin", "X:=", str(origin[0]), "Y:=", str(origin[1]), "Z:=", str(origin[2])],
                                ["NAME:X Axis", "X:=", str(new_x[0]), "Y:=", str(new_x[1]), "Z:=", str(new_x[2])],
                                ["NAME:Y Point", "X:=", str(new_y[0]), "Y:=", str(new_y[1]), "Z:=", str(new_y[2])]]]])

    def set_coor_sys(self, coor_sys):
        if coor_sys != self.get_coor_sys():
            self._modeler.SetWCS(["NAME:SetWCS Parameter","Working Coordinate System:=", coor_sys , "RegionDepCSOk:="	, False])

    def get_coor_sys(self):
        return self._modeler.GetActiveCoordinateSystem()

    def delete(self, entity):
        objects = [self._modeler.GetObjectName(str(ii))
                        for ii in range(int(self._modeler.GetNumObjects()))]
        if entity.name in objects:
            self._modeler.Delete(["NAME:Selections",
                                  "Selections:=", entity.name])

    def delete_all_objects(self):
        objects = [self._modeler.GetObjectName(str(ii))
                        for ii in range(int(self._modeler.GetNumObjects()))]
        self._modeler.Delete(self._selections_array(*objects))

    @assert_name
    def box(self, pos, size, **kwargs):
        if len(pos)==2:
            pos.append(0)
        if len(size)==2:
            size.append(0)
        pos = parse_entry(pos)
        size = parse_entry(size)
        name = self._modeler.CreateBox(
            ["NAME:BoxParameters",
             "XPosition:=", str(pos[0]),
             "YPosition:=", str(pos[1]),
             "ZPosition:=", str(pos[2]),
             "XSize:=", str(size[0]),
             "YSize:=", str(size[1]),
             "ZSize:=", str(size[2])],
            self._attributes_array(**kwargs))
        return name

    @assert_name
    def box_center(self, pos, size, **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        corner_pos = [val(p) - val(s)/2 for p, s in zip(pos, size)]
        return self.box(corner_pos, size, **kwargs)

    @assert_name
    def polyline(self, points, closed=True, **kwargs):
        for i in range(len(points)):
            if isinstance(points[i], tuple) and len(points[i])==2:
                points[i]+=(0,)
            elif isinstance(points[i], list) and len(points[i])==2:
                points[i].append(0)

        points = parse_entry(points)
        pointsStr = ["NAME:PolylinePoints"]
        indexsStr = ["NAME:PolylineSegments"]
#        print("final_points", points)
        for ii, point in enumerate(points):
            pointsStr.append(["NAME:PLPoint", "X:=", str(point[0]), "Y:=", str(point[1]), "Z:=", str(point[2])])
            indexsStr.append(["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", ii, "NoOfPoints:=", 2])
        if closed:
            pointsStr.append(["NAME:PLPoint", "X:=", str(points[0][0]), "Y:=", str(points[0][1]), "Z:=", str(points[0][2])])
            params_closed = ["IsPolylineCovered:=", True, "IsPolylineClosed:=", True]
        else:
            indexsStr = indexsStr[:-1]
            params_closed = ["IsPolylineCovered:=", True, "IsPolylineClosed:=", False]


        name = self._modeler.CreatePolyline(
            ["NAME:PolylineParameters",
            *params_closed,
            pointsStr,
            indexsStr],
            self._attributes_array(**kwargs))

        return name

    @assert_name
    def rect(self, pos, size, **kwargs):
        if len(pos)==2:
            pos.append(0)
        if len(size)==2:
            size.append(0)
        pos = parse_entry(pos)
        size = parse_entry(size)
        index = size.index(0)
        if index>=0:
            axis = "XYZ"[index]
        w_idx, h_idx = {'X': (1, 2),
                        'Y': (2, 0),
                        'Z': (0, 1)}[axis]
        name = self._modeler.CreateRectangle(
            ["NAME:RectangleParameters",
             "XStart:=", str(pos[0]),
             "YStart:=", str(pos[1]),
             "ZStart:=", str(pos[2]),
             "Width:=", str(size[w_idx]),
             "Height:=", str(size[h_idx]),
             "WhichAxis:=", axis],
            self._attributes_array(**kwargs)
        )
        return name

    @assert_name
    def rect_center(self, pos, size, **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        corner_pos = [val(p) - val(s)/2 for p, s in zip(pos, size)]
        name = self.rect(corner_pos, size, **kwargs)
        return name

    @assert_name
    def sphere(self, pos, radius, **kwargs):
        name = self._modeler.CreateSphere(
            ["NAME:SphereParameters",
             "XCenter:=", str(pos[0]),
             "YCenter:=", str(pos[1]),
             "ZCenter:=", str(pos[2]),
             "Radius:=", str(radius)],
            self._attributes_array(**kwargs))
        return name
    
    @assert_name
    def torus(self, pos, majorradius, minorradius, axis, **kwargs):
        assert axis in "XYZ"
        name = self._modeler.CreateTorus(
            ["NAME:TorusParameters",
             "XCenter:=", str(pos[0]),
             "YCenter:=", str(pos[1]),
             "ZCenter:=", str(pos[2]),
             "MajorRadius:=", str(majorradius),
             "MinorRadius:=", str(minorradius),
             "WhichAxis:=", axis],
            self._attributes_array(**kwargs))
        return name
        
    @assert_name
    def cylinder(self, pos, radius, height, axis, **kwargs):
        assert axis in "XYZ"
        name = self._modeler.CreateCylinder(
            ["NAME:CylinderParameters",
             "XCenter:=", str(pos[0]),
             "YCenter:=", str(pos[1]),
             "ZCenter:=", str(pos[2]),
             "Radius:=", str(radius),
             "Height:=", str(height),
             "WhichAxis:=", axis,
             "NumSides:=", 0],
            self._attributes_array(**kwargs))
        return name

    @assert_name
    def cylinder_center(self, pos, radius, height, axis, **kwargs):
        assert axis in "XYZ"
        axis_idx = ["X", "Y", "Z"].index(axis)
        edge_pos = copy(pos)
        edge_pos[axis_idx] = val(pos[axis_idx]) - val(height)/2
        return self.cylinder(edge_pos, radius, height, axis, **kwargs)

    @assert_name
    def disk(self, pos, radius, axis, **kwargs):
        assert axis in "XYZ"
        name = self._modeler.CreateEllipse(
            ["NAME:EllipsdeParameters",
             "XCenter:=", str(pos[0]),
             "YCenter:=", str(pos[1]),
             "ZCenter:=", str(pos[2]),
             "MajRadius:=", str(radius),
             "Ratio:=", 1,
             "WhichAxis:=", axis],
             self._attributes_array(**kwargs))
        return name

    @assert_name
    def wirebond(self, pos, ori, ymax, ymin, height='0.1mm', **kwargs): #ori should be normed
        bond_diam = '20um'
        params = parse_entry(pos, ori, ymax, ymin, height, bond_diam)
        pos, ori, ymax, ymin, heigth, bond_diam = params
        bond1 = pos + ori.orth()*(ymax+2*bond_diam)
        width = ymax-ymin + 4*bond_diam
        direction = - ori.orth()
        xpad, ypad, _ = bond1
        xdir, ydir, _ = direction
        kwargs['material'] = 'perfect conductor'
        kwargs['solve_inside'] = False
        name = self._modeler.CreateBondwire(["NAME:BondwireParameters",
                                            "WireType:=", "Low",
                                            "WireDiameter:=", str(bond_diam),
                                            "NumSides:=", str(6),
                                            "XPadPos:=", str(xpad),
                                            "YPadPos:=", str(ypad),
                                            "ZPadPos:=", "0mm",
                                            "XDir:=", str(xdir),
                                            "YDir:=", str(ydir),
                                            "ZDir:=", str(0),
                                            "Distance:=", str(width),
                                            "h1:=", str(height),
                                            "h2:=", "0mm",
                                            "alpha:=", "80deg",
                                            "beta:=", "80deg",
                                            "WhichAxis:=", "Z"],
                                            self._attributes_array(**kwargs))
        return name

    def sweep_along_path(self, entity_to_sweep, path_entity):
#        name_temp = path_entity.name
#        self.rename_entity(path_entity, path_entity.name+'_path')
#        self.rename_entity(entity_to_sweep, name_temp)

        names = [entity_to_sweep.name, path_entity.name]
        self._modeler.SweepAlongPath(self._selections_array(*names),
                                     ["NAME:PathSweepParameters",
                                		"DraftAngle:="		, "0deg",
                                		"DraftType:="		, "Round",
                                		"CheckFaceFaceIntersection:=", False,
                                		"TwistAngle:="		, "0deg"])
        entity_to_sweep.dimension += 1
        return entity_to_sweep.name

    def duplicate_along_line(self, entity, vec, n=2):
        self._modeler.DuplicateAlongLine(["NAME:Selections","Selections:=", entity.name,
                                          "NewPartsModelFlag:="	, "Model"],
                                        	["NAME:DuplicateToAlongLineParameters",
                                    		"CreateNewObjects:="	, True,
                                    		"XComponent:="		, str(vec[0]),
                                    		"YComponent:="		, str(vec[1]),
                                    		"ZComponent:="		, str(vec[2]),
                                    		"NumClones:="		, str(n)],
                                        	["NAME:Options",
                                        	"DuplicateAssignments:=", False],
                                        	["CreateGroupsForNewObjects:=", False	])


    def intersect(self, entities, keep_originals=False):
        names = []
        for entity in entities:
            if entity!=None:
                names.append(entity.name)
        self._modeler.Intersect(
                                self._selections_array(*names),
                                ["NAME:IntersectParameters", "KeepOriginals:=", keep_originals]
                                )
        return names[0]

#    def separate_bodies(self, name):
#        self._modeler.SeparateBody(["NAME:Selections",
#                                		"Selections:=", name,
#                                		"NewPartsModelFlag:="	, "Model"
#                                	],
#                                	[
#                                		"CreateGroupsForNewObjects:=", False
#                                	])

    def assign_perfect_E(self, entities, name):
        if not isinstance(entities, list):
            entities = [entities]
        entity_names = [entity.name for entity in entities]
        self._boundaries.AssignPerfectE(["NAME:"+name,
                                         "Objects:=", entity_names,
                                         "InfGroundPlane:=", False])

    def assign_perfect_E_faces(self, entity):
        # this is very peculiar to cavity Si chips
        faces = list(self.get_face_ids(entity.name))
        faces = [int(ii) for ii in faces]
        faces.sort()
        faces_perfE = [faces[1]]+faces[6:]
        faces_loss = faces[2:6]

        self._boundaries.AssignPerfectE(
        	[
        		"NAME:"+str(name)+'_PerfE',
        		"Faces:="		, faces_perfE,
        		"InfGroundPlane:="	, False
        	])

        self._boundaries.AssignFiniteCond([ "NAME:FiniteCond3",
                                        		"Faces:="		, faces_loss,
                                        		"UseMaterial:="		, True,
                                        		"Material:="		, "lossy conductor",
                                        		"UseThickness:="	, False,
                                        		"Roughness:="		, "0um",
                                        		"InfGroundPlane:="	, False,
                                        		"IsTwoSided:="		, False,
                                        		"IsInternal:="		, True ])

    def assign_waveport(self, entity, name):
        Nmodes=1
        PolarizeEField=False
        RenormValue="50Ohm"
        UseIntLine=False
        DoDeembed=False
        # creates a waveport on a face
        faces = list(self.get_face_ids(entity.name))
        faces = [int(ii) for ii in faces]
        # creates name for modes
        modesarray = ["NAME:Modes"]
        for n in range(Nmodes):
            modesarray.append(["NAME:Mode" + str(n + 1),
                               "ModeNum:=", n + 1,
                               "UseIntLine:=", UseIntLine])

        self._boundaries.AssignWavePort(["NAME:", name,
                                         "NumModes:=", Nmodes,
                                         "PolarizeEField:=", PolarizeEField,
                                         "DoDeembed:=", DoDeembed,
                                         "RenormValue:=", RenormValue, 
                                          modesarray, 
                                         "Faces:=", faces])

    def assign_mesh_length(self, entities, length):
        if not isinstance(entities, list):
            entities = [entities]
        name = entities[0].name
        params = ["NAME:"+name]
        list_bool = [entity.dimension==3 for entity in entities]
        if all(list_bool):
            params += ["RefineInside:=", True, "Enabled:=", True]
        elif any(list_bool):
            raise TypeError('Please do not mix 3D element with other \
                            dimenstions when assigning mesh length')
        else:
            params += ["RefineInside:=", False, "Enabled:=", True]

        params += ["Objects:=", [entity.name for entity in entities]]
        params += ["RestrictElem:=", False,
			         "RestrictLength:=",  True,
			         "MaxLength:=", str(length)]
        self._mesh.AssignLengthOp(params)

    def assign_lumped_rlc(self, entity, r, l, c, start, end, name="RLC"):
        # parsing start and end is very peculiar for this function
        # hfss style
        start, end = parse_entry(start, end)
        start, end = val(start, end)
        start = [str(elt*1000)+'mm' for elt in start]
        end = [str(elt*1000)+'mm' for elt in end]
        params = ["NAME:"+entity.name+'_'+name, "Objects:="]
        params.append([entity.name])
        params.append(["NAME:CurrentLine", "Start:=", start, "End:=", end])
        params += ["UseResist:=", r != 0, "Resistance:=", str(r),
                   "UseInduct:=", l != 0, "Inductance:=", str(l),
                   "UseCap:=", c != 0, "Capacitance:=", str(c)]
        self._boundaries.AssignLumpedRLC(params)


    def create_object_from_face(self, entity):
        faces = list(self._modeler.GetFaceIDs(entity.name))
        faces.sort()
        face = faces[0]
        self._modeler.CreateObjectFromFaces(["NAME:Selections",
                                            "Selections:="		, entity.name,
                                            "NewPartsModelFlag:="	, "Model"],
            ["NAME:Parameters",["NAME:BodyFromFaceToParameters","FacesToDetach:="	, [int(face)]]],
            ["CreateGroupsForNewObjects:=", False])
        return entity.name+'_ObjectFromFace1'

    def eval_expr(self, expr, units="mm"):
        if not isinstance(expr, str):
            return expr
        return self.parent.eval_expr(expr, units)

    def eval_var_str(self, name, unit=None):
        raise NotImplementedError('eval_var_str is not implemented anymore')
        # if not isinstance(name, VariableString):
        #     if unit is not None:
        #         return str(name)+' '+unit
        #     else:
        #         return str(name)
        # return self.parent.eval_var_str(name, unit=unit)

    def fillet(self, entity, radius, vertex_indices=None):
        # if the geometry is simple : rect, polyline, fillet indexing is consistent
        # if the geometry was subtracted/united : fillet indexing is not consistent
        # if the geometry ha a hole, should not use fillet, since un expected behaviour
        # if object has already been filleted should not be possible to fillet twice
        hfss_vertex_indices = self.get_vertex_ids(entity)
        if vertex_indices is None:
            if entity.dimension == 1:
                hfss_vertex_indices = hfss_vertex_indices[1:-1]
            to_fillet = [int(vertex) for vertex in hfss_vertex_indices]
        else:
            if not isinstance(vertex_indices, list):
                vertex_indices = [vertex_indices]
            to_fillet = [int(hfss_vertex_indices[v]) for v in vertex_indices]
        if to_fillet != []:
            self._modeler.Fillet(["NAME:Selections", "Selections:=", entity.name],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", [],
                                "Vertices:=", to_fillet,
                                "Radius:=", str(radius),
                                "Setback:=", "0mm"]])

    def _fillet_edges(self, entity, radius, edge_index):
        edges = self.get_edge_ids(entity.name)
        if isinstance(edge_index, list):
            to_fillet = [int(edges[e]) for e in edge_index]
        else:
            to_fillet = [int(edges[edge_index])]
        self._modeler.Fillet(["NAME:Selections", "Selections:=", entity.name],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", to_fillet,
                                "Vertices:=", [],
                                "Radius:=", str(radius),
                                "Setback:=", "0mm"]])

    def get_faces(self, entity):
        '''
        Only for boxes
        Output :
           x_back_face,
           x_front_face,
           y_back_face,
           y_front_face,
           z_back_face,
           z_front_face
        '''
        assert entity.dimension==3, 'This function is only callable for 3D objects'
        faces = self.get_face_ids(entity)
        z_back_face, z_front_face = faces[0], faces[1]
        y_back_face, y_front_face = faces[2], faces[4]
        x_back_face, x_front_face = faces[3], faces[5]
        return x_back_face, x_front_face, y_back_face, y_front_face, z_back_face, z_front_face

    def get_face_ids(self, entity):
        return self._modeler.GetFaceIDs(entity.name)

    def get_vertex_ids(self, entity):
        return self._modeler.GetVertexIDsFromObject(entity.name)

    def get_vertices(self, entity):
        vertices_ids = self.get_vertex_ids(entity)
        return [[*map(float, self._modeler.GetVertexPosition(vertex)[:2])]
                for vertex in vertices_ids]

    def get_edge_ids(self, entity):
        return self._modeler.GetEdgeIDsFromObject(entity.name)

    def get_matched_object_name(self, name):
        return self._modeler.GetMatchedObjectName(name+'*')


    def mirrorZ(self, entity):
        self._modeler.Mirror([
                        		"NAME:Selections", "Selections:=", entity.name,
                        		"NewPartsModelFlag:="	, "Model"
                        	  ],
                        	  [
                        		"NAME:MirrorParameters",
                        		"MirrorBaseX:="		, "0mm",
                        		"MirrorBaseY:="		, "0mm",
                        		"MirrorBaseZ:="		, "0mm",
                        		"MirrorNormalX:="	, "0mm",
                        		"MirrorNormalY:="	, "0mm",
                        		"MirrorNormalZ:="	, "1mm"
                        	  ])
        return entity.name

    def make_center_line(self, entity, axis):
        raise NotImplementedError()
        center = [c + s/2 if s else c for c, s in zip(entity.pos, entity.size)]
        axis_idx = ["x", "y", "z"].index(axis.lower())
        start = [c for c in center]
        start[axis_idx] -= entity.size[axis_idx]/2
        start = [self.eval_var_str(s, unit=LENGTH_UNIT) for s in start]
        end = [c for c in center]
        end[axis_idx] += entity.size[axis_idx]/2
        end = [self.eval_var_str(s, unit=LENGTH_UNIT) for s in end]
        return start, end

    def make_rlc_boundary(self, entity, axis, r=0, l=0, c=0, name="LumpRLC"):
        name = entity.name+'_'+name
        start, end = self.make_center_line(entity, axis)
        self.modeler._make_lumped_rlc(entity, r, l, c, start, end, name=name)

    def _make_lumped_rlc(self, entity, r, l, c, start, end, name="LumpLRC"):
        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:"+name, "Objects:="]
        params.append([entity.name])
        params.append(["NAME:CurrentLine", "Start:=", start, "End:=", end])
        params += ["UseResist:=", r != 0, "Resistance:=", r,
                   "UseInduct:=", l != 0, "Inductance:=", l,
                   "UseCap:=", c != 0, "Capacitance:=", c]
        self._boundaries.AssignLumpedRLC(params)

    def _make_lumped_port(self, start, end, obj_arr, z0="50ohm", name="LumpPort"):
        #TODO
        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:"+name]
        params += obj_arr
        params += ["RenormalizeAllTerminals:=", True, "DoDeembed:=", False,
                   ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,
                                   ["NAME:IntLine", "Start:=", start, "End:=", end],
                                   "CharImp:=", "Zpi", "AlignmentGroup:=", 0, "RenormImp:=", "50ohm"]],
                   "ShowReporterFilter:=", False, "ReporterFilter:=", [True],
                   "FullResistance:=", "50ohm", "FullReactance:=", "0ohm"]

        self._boundaries.AssignLumpedPort(params)

    def assign_material(self, entities, material):
        if not isinstance(entities, list):
            entities = [entities]
        entities_name = [entity.name for entity in entities]
        solve_inside = not(material=='perfect conductor')
        material = '\"'+material+'\"'
        self._modeler.AssignMaterial(self._selections_array(*entities_name),
                        	[
                        		"NAME:Attributes",
                        		"MaterialValue:="	, material,
                        		"SolveInside:="		, solve_inside,
                        		"IsMaterialEditable:="	, True,
                        		"UseMaterialAppearance:=", False
                        	])

    def rotate(self, entities, angle):
        if not isinstance(entities, list):
            entities = [entities]
        names = [entity.name for entity in entities]
        self._modeler.Rotate(self._selections_array(*names),
            ["NAME:RotateParameters", "RotateAxis:=", "Z", "RotateAngle:=",
             "%ddeg"%(angle)])

    def rename(self, entity, name):
        new_name = self._modeler.ChangeProperty(["NAME:AllTabs",
      		["NAME:Geometry3DAttributeTab",
      		["NAME:PropServers", str(entity.name)],
      			["NAME:ChangedProps",["NAME:Name","Value:=", str(name)]]]])
        return new_name

    def set_units(self, units='m'):
        self._modeler.SetModelUnits(["NAME:Units Parameter",
                                     "Units:=", units,"Rescale:=",False])

    def subtract(self, blank_entities, tool_entities, keep_originals=False):
        blank_names = []
        for entity in blank_entities:
            blank_names.append(entity.name)
        tool_names = []
        for entity in tool_entities:
            tool_names.append(entity.name)
        selection_array= ["NAME:Selections",
                          "Blank Parts:=", ",".join(blank_names),
                          "Tool Parts:=", ",".join(tool_names)]
        self._modeler.Subtract(selection_array,
                                ["NAME:UniteParameters",
                                 "KeepOriginals:=", keep_originals])



    def sweep_along_vector(self, entities, vector):
        names = [entity.name for entity in entities]
        self._modeler.SweepAlongVector(self._selections_array(*names),
                                        	["NAME:VectorSweepParameters",
                                        		"DraftAngle:="		, "0deg",
                                        		"DraftType:="		, "Round",
                                        		"CheckFaceFaceIntersection:=", False,
                                        		"SweepVectorX:="	, str(vector[0]),
                                        		"SweepVectorY:="	, str(vector[1]),
                                        		"SweepVectorZ:="	, str(vector[2])
                                        	])

    #typically use for trenching of the gapobjects.
    def thicken_sheet(self, entity, thickness, bothsides=False):
        self._modeler.ThickenSheet([
                                		"NAME:Selections", "Selections:=", entity.name,
                                		"NewPartsModelFlag:="	, "Model"
                                	],
                                	[
                                		"NAME:SheetThickenParameters",
                                		"Thickness:="		, str(thickness),
                                		"BothSides:="		, bothsides
                                	])
    def translate(self, entities, vector):
        if not isinstance(entities, list):
            entities = [entities]
        names = [entity.name for entity in entities]
        self._modeler.Move(self._selections_array(*names),
                            ["NAME:TranslateParameters",
                        		"TranslateVectorX:="	, str(vector[0]),
                        		"TranslateVectorY:="	, str(vector[1]),
                        		"TranslateVectorZ:="	, str(vector[2])])

    def unite(self, entities, keep_originals=False):
        names = [entity.name for entity in entities]
        self._modeler.Unite(self._selections_array(*names),
            ["NAME:UniteParameters", "KeepOriginals:=", keep_originals])
        return entities.pop(0)



# class ModelEntity():
#     instances_layered = {}
#     dict_instances = {}
#     instances_to_move = []
#     def __init__(self, name, dimension, coor_sys, model = 'True', layer=layer_Default):# model,
#         self.name = name
#         self.dimension = dimension
#         self.coor_sys = coor_sys
#         self.history = []
#         self.model = model
#         self.layer = layer

#         ModelEntity.dict_instances[name] = self
#         ModelEntity.find_last_list(ModelEntity.instances_to_move).append(self)
#         if layer in self.instances_layered.keys():
#             ModelEntity.instances_layered[layer].append(self)
#         else:
#             ModelEntity.instances_layered[layer]=[self]

#     def __str__(self):
#         return self.name

#     def check_name(self, name):
#         i = 0
#         new_name = name
#         while(new_name in self.dict_instances.keys()):
#             i+=1
#             new_name = name+'_'+str(i)
#         return new_name

#     @classmethod
#     def printInstances(cls):
#         for instance in cls.instances:
#             print(instance)

#     @staticmethod
#     def find_last_list(list_entities=instances_to_move):
# #        print(list_entities)
# #        if isinstance(list_entities, Port) :
# #            return None
#         if isinstance(list_entities, list):
#             if len(list_entities)==0:
#                 return list_entities
#             else:
#                 if isinstance(list_entities[-1], list):
#                     return ModelEntity.find_last_list(list_entities[-1])
#                 else:
#                     return list_entities
#         else:
#             return list_entities


# #    @staticmethod
# #    def merge(n):
# #        key_list = []
# #        for idx, key in enumerate(ModelEntity.instances_to_move):
# #            if idx>=n:
# #                key_list.append(key)
# #        for key in key_list:
# #            entity = ModelEntity.instances_to_move.pop(key, None)
# #            ModelEntity.instances_moved[key]=entity

#     @staticmethod
#     def reset():
#         ModelEntity.instances_layered = {}
#         ModelEntity.dict_instances = {}
#         ModelEntity.instances_to_move = []

#     def append_history(self, entity):
#         self.history.append(entity)

#     def delete(self):
#         del self

#     def copy(self, dimension=None):
#         #TODO complete
#         if dimension==None:
#             dimension = self.dimension
#         return ModelEntity(self.name, dimension, self.coor_sys)

#     def rename_entity(self, name):
#         self.name = name

#     def modify_dimension(self, dimension):
#         self.dimension = dimension

#     def _sweep_along_path(self, path_obj):
#         self.dimension = 2
#         self.append_history(path_obj)
#         return self

#     def sweep_along_vector(self, vect):
#         "Duplicate the line and create a surface with the copy"
#         self.dimension = 2
#         return self

#     def thicken_sheet(self, thickness, bothsides=False):
#         self.dimension = 3
#         return self

class HfssFieldsCalc(COMWrapper):
    def __init__(self, setup):
        """
        :type setup: HfssSetup
        """
        super(HfssFieldsCalc, self).__init__()
        self.parent = setup
        self.Mag_E = NamedCalcObject("Mag_E", setup)
        self.Mag_H = NamedCalcObject("Mag_H", setup)
        self.Mag_Jsurf = NamedCalcObject("Mag_Jsurf", setup)
        self.Mag_Jvol = NamedCalcObject("Mag_Jvol", setup)
        self.Vector_E = NamedCalcObject("Vector_E", setup)
        self.Vector_H = NamedCalcObject("Vector_H", setup)
        self.Vector_Jsurf = NamedCalcObject("Vector_Jsurf", setup)
        self.Vector_Jvol = NamedCalcObject("Vector_Jvol", setup)
        self.ComplexMag_E = NamedCalcObject("ComplexMag_E", setup)
        self.ComplexMag_H = NamedCalcObject("ComplexMag_H", setup)
        self.ComplexMag_Jsurf = NamedCalcObject("ComplexMag_Jsurf", setup)
        self.ComplexMag_Jvol = NamedCalcObject("ComplexMag_Jvol", setup)
        self.P_J = NamedCalcObject("P_J", setup)

    def clear_named_expressions(self):
        self.parent.parent._fields_calc.ClearAllNamedExpr()

class CalcObject(COMWrapper):
    def __init__(self, stack, setup):
        """
        :type stack: [(str, str)]
        :type setup: HfssSetup
        """
        super(CalcObject, self).__init__()
        self.stack = stack
        self.setup = setup
        self.calc_module = setup.parent._fields_calc

    def _bin_op(self, other, op):
        if isinstance(other, (int, float)):
            other = ConstantCalcObject(other, self.setup)

        stack = self.stack + other.stack
        stack.append(("CalcOp", op))
        return CalcObject(stack, self.setup)

    def _unary_op(self, op):
        stack = self.stack[:]
        stack.append(("CalcOp", op))
        return CalcObject(stack, self.setup)

    def __add__(self, other):
        return self._bin_op(other, "+")

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self._bin_op(other, "-")

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        return self._bin_op(other, "*")

    def __rmul__(self, other):
        return self*other

    def __div__(self, other):
        return self._bin_op(other, "/")

    def __rdiv__(self, other):
        other = ConstantCalcObject(other, self.setup)
        return other/self

    def __pow__(self, other):
        return self._bin_op(other, "Pow")

    def dot(self, other):
        return self._bin_op(other,"Dot")

    def __neg__(self):
        return self._unary_op("Neg")

    def __abs__(self):
        return self._unary_op("Abs")

    def __mag__(self):
        return self._unary_op("Mag")

    def conj(self):
        return self._unary_op("Conj") # make this right

    def scalar_x(self):
        return self._unary_op("ScalarX")

    def scalar_y(self):
        return self._unary_op("ScalarY")

    def scalar_z(self):
        return self._unary_op("ScalarZ")

    def norm_2(self):
        return (self.__mag__()).__pow__(2)
        #return self._unary_op("ScalarX")**2+self._unary_op("ScalarY")**2+self._unary_op("ScalarZ")**2

    def real(self):
        return self._unary_op("Real")

    def imag(self):
        return self._unary_op("Imag")

    def _integrate(self, name, type):
        stack = self.stack + [(type, name), ("CalcOp", "Integrate")]
        return CalcObject(stack, self.setup)

    def getQty(self, name):
        stack = self.stack + [("EnterQty", name)]
        return CalcObject(stack, self.setup)

    def integrate_line(self, name):
        return self._integrate(name, "EnterLine")

    def integrate_line_tangent(self, name):
        ''' integrate line tangent to vector expression \n
            name = of line to integrate over '''
        self.stack = self.stack + [("EnterLine", name),
                                   ("CalcOp",    "Tangent"),
                                   ("CalcOp",    "Dot")]#,
                              #("EnterLine", name),
                              #("CalcOp",   "Integrate")]
        return self.integrate_line(name)

    def integrate_surf(self, name="AllObjects"):
        return self._integrate(name, "EnterSurf")

    def integrate_vol(self, name="AllObjects"):
        return self._integrate(name, "EnterVol")

    def times_eps(self):
        stack = self.stack + [("ClcMaterial", ("Permittivity (epsi)", "mult"))]
        return CalcObject(stack, self.setup)

    def times_mu(self):
        stack = self.stack + [("ClcMaterial", ("Permeability (mu)", "mult"))]
        return CalcObject(stack, self.setup)

    def write_stack(self):
        for fn, arg in self.stack:
            if numpy.size(arg)>1:
                getattr(self.calc_module, fn)(*arg)
            else:
                getattr(self.calc_module, fn)(arg)

    def save_as(self, name):
        """if the object already exists, try clearing your
        named expressions first with fields.clear_named_expressions"""
        self.write_stack()
        self.calc_module.AddNamedExpr(name)
        return NamedCalcObject(name, self.setup)

    def evaluate(self, phase=0, lv=None, print_debug = False):#, n_mode=1):
        self.write_stack()
        if print_debug:
            print('---------------------')
            print('writing to stack: OK')
            print('-----------------')
        #self.calc_module.set_mode(n_mode, 0)
        setup_name = self.setup.solution_name

        if lv is not None:
           args = lv
        else:
           args = []

        args.append("Phase:=")
        args.append(str(int(phase)) + "deg")

        if isinstance(self.setup, HfssDMSetup):
            args.extend(["Freq:=", self.setup.solution_freq])

        self.calc_module.ClcEval(setup_name, args)
        return float(self.calc_module.GetTopEntryValue(setup_name, args)[0])

class NamedCalcObject(CalcObject):
    def __init__(self, name, setup):
        self.name = name
        stack = [("CopyNamedExprToStack", name)]
        super(NamedCalcObject, self).__init__(stack, setup)

class ConstantCalcObject(CalcObject):
    def __init__(self, num, setup):
        stack = [("EnterScalar", num)]
        super(ConstantCalcObject, self).__init__(stack, setup)

def get_desktop():
    import ctypes, os
    try:
        is_admin = os.getuid() == 0
    except AttributeError:
        is_admin = ctypes.windll.shell32.IsUserAnAdmin() != 0
    if not is_admin:
        print('\033[93m WARNING: you are not runnning as an admin! You need to run as an admin. You will probably get an error next. \033[0m')

    app = HfssApp()
    return app.get_app_desktop()

def get_active_project():
    ''' If you see the error:
        "The requested operation requires elevation."
        then you need to run your python as an admin.
    '''
    import ctypes, os
    try:
        is_admin = os.getuid() == 0
    except AttributeError:
        is_admin = ctypes.windll.shell32.IsUserAnAdmin() != 0
    if not is_admin:
        print('\033[93m WARNING: you are not runnning as an admin! You need to run as an admin. You will probably get an error next. \033[0m')

    app = HfssApp()
    desktop = app.get_app_desktop()
    return desktop.get_active_project()

def get_active_design():
    project = get_active_project()
    return project.get_active_design()

def get_report_arrays(name):
    d = get_active_design()
    r = HfssReport(d, name)
    return r.get_arrays()

def load_HFSS_project(proj_name, project_path, extension = '.aedt'):  # aedt is for 2016 version
    ''' proj_name == None => get active.
        (make sure 2 run as admin)
    '''

    # Checks
    assert os.path.isdir(project_path),   "ERROR! project_path is not a valid directory. Check the path, and especially \\ charecters."

    project_path +=  proj_name + extension

    if os.path.isfile(project_path):
        print('\tFile path to HFSS project found.')
    else:
        raise Exception("ERROR! Valid directory, but invalid project filename. Not found! Please check your filename.\n%s\n" % project_path )

    if os.path.isfile(project_path+'.lock'):
        print('\t\tFile is locked. If connection fails, delete the .lock file.'    )

    app     = HfssApp()
    print("\tOpened HfssApp.")
    desktop = app.get_app_desktop()
    print( "\tOpened Hfss desktop.")
    if proj_name is not None:
        if proj_name in desktop.get_project_names():
            desktop.set_active_project(proj_name)
            project = desktop.get_active_project()
        else:
            project = desktop.open_project(project_path)
    else:
        project = desktop.get_active_project()
    print("\tOpened link to HFSS user project.")
    return app, desktop, project

