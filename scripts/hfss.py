from __future__ import division
import atexit
from copy import copy
import os

import dis
import inspect
import tempfile
import types
import numpy
import signal
import pythoncom
import time

from sympy.parsing import sympy_parser
from pint import UnitRegistry
from win32com.client import Dispatch, CDispatch

ureg = UnitRegistry()
Q = ureg.Quantity

BASIS_ORDER = {"Zero Order": 0,
               "First Order": 1,
               "Second Order": 2,
               "Mixed Order": -1}
LENGTH = '[length]'
INDUCTANCE = '[length] ** 2 * [mass] / [current] ** 2 / [time] ** 2'
CAPACITANCE = '[current] ** 2 * [time] ** 4 / [length] ** 2 / [mass]'
RESISTANCE = '[length] ** 2 * [mass] / [current] ** 2 / [time] ** 3'

LENGTH_UNIT = 'meter'
INDUCTANCE_UNIT = 'nH'
CAPACITANCE_UNIT = 'fF'
RESISTANCE_UNIT = 'ohm'


def simplify_arith_expr(expr):
    try:
        out = repr(sympy_parser.parse_expr(str(expr)))
        return out
    except Exception:
        print("Couldn't parse", expr)
        raise


def increment_name(base, existing):
    if base not in existing:
        return base
    n = 1
    make_name = lambda: base + str(n)
    while make_name() in existing:
        n += 1
    return make_name()


def extract_value_unit(expr, units):
    """
    :type expr: str
    :type units: str
    :return: float
    """
    try:
        return Q(expr).to(units).magnitude
    except Exception:
        try:
            return float(expr)
        except Exception:
            return expr


def extract_value_dim(expr):
    """
    type expr: str
    """
    return str(Q(expr).dimensionality)


def parse_entry(entry):
    #should take a list of tuple of list... of int, float or str...
    if not isinstance(entry, list) and not isinstance(entry, tuple):
        return extract_value_unit(entry, LENGTH_UNIT)
    else:
        entries = entry
        _entry = []
        for entry in entries:
            _entry.append(parse_entry(entry))
        return _entry


class VariableString(str):
    # TODO: What happen with a list (Vector in our case)
    def __add__(self, other):
        other = self.rem_unit(other)
        if other=="'":
            return super(VariableString, self, other).__add__()
        return var("(%s) + (%s)" % (self, other))

    def __radd__(self, other):
        other = self.rem_unit(other)
        if other=="'":
            return super(VariableString, self, other).__radd__()
        return var("(%s) + (%s)" % (other, self))

    def __sub__(self, other):
        other = self.rem_unit(other)
        return var("(%s) - (%s)" % (self, other))

    def __rsub__(self, other):
        other = self.rem_unit(other)
        return var("(%s) - (%s)" % (other, self))

    def __mul__(self, other):
        other = self.rem_unit(other)
        return var("(%s) * (%s)" % (self, other))

    def __rmul__(self, other):
        other = self.rem_unit(other)
        return var("(%s) * (%s)" % (other, self))

    def __div__(self, other):
        other = self.rem_unit(other)
        return var("(%s) / (%s)" % (self, other))

    def __rdiv__(self, other):
        other = self.rem_unit(other)
        return var("(%s) / (%s)" % (other, self))

    def __truediv__(self, other):
        other = self.rem_unit(other)
        return var("(%s) / (%s)" % (self, other))

    def __rtruediv__(self, other):
        other = self.rem_unit(other)
        return var("(%s) / (%s)" % (other, self))

    def __pow__(self, other):
        other = self.rem_unit(other)
        return var("(%s) ** (%s)" % (self, other))

#    def __rpow__(self, other):
#        other = self.rem_unit(other)
#        return var("(%s) ** (%s)" % (other, self))

    def __neg__(self):
        return var("-(%s)" % self)

    def __abs__(self):
        return var("abs(%s)" % self)

    def rem_unit(self, other):
        try:
            _other = extract_value_unit(other, LENGTH_UNIT)
            return _other
        except Exception:
            return other


def var(x):
    if isinstance(x, str):
        return VariableString(simplify_arith_expr(x))
    return x


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
        return [VariableString(s) for s in self._project.GetVariables()]

    def get_variables(self):
        return {VariableString(s): self.get_variable_value(s) for s in self._project.GetVariables()}

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
        return VariableString(name)

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
                "Value:=", value]]]])
        self.store_variable(name, value)
        return VariableString(name)

    def set_variable(self, name, value, postprocessing=False):
        # TODO: check if variable does not exist and quit if it doesn't?
#        print(name)
#        print(type(name))
#        print(value)
#        print(type(value))
        if name not in self._design.GetVariables()+self._design.GetPostProcessingVariables():
            self.create_variable(name, value, postprocessing=postprocessing)
        else:
            self._design.SetVariableValue(name, value)
            self.store_variable(name, value)
        return VariableString(name)

    def store_variable(self, name, value):
        if not isinstance(value, VariableString):
            if LENGTH == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, LENGTH_UNIT)
            if INDUCTANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, INDUCTANCE_UNIT)
            if CAPACITANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, CAPACITANCE_UNIT)
            if RESISTANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, RESISTANCE_UNIT)
        else:
            self.variables[name]=value

    def get_variable_value(self, name):
        return self._design.GetVariableValue(name)

    def get_variable_names(self):
        return [VariableString(s) for s in self._design.GetVariables()+self._design.GetPostProcessingVariables()]

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

    def set_units(self, units='m'):
        self._modeler.SetModelUnits(["NAME:Units Parameter",
                                     "Units:=", units,"Rescale:=",False])

    def _attributes_array(self, name=None, nonmodel=False, color=None,
                          transparency=0.9, material=None, solve_inside = None):
        arr = ["NAME:Attributes", "PartCoordinateSystem:=", "Global"]
        if name is not None:
            arr.extend(["Name:=", name])
        if nonmodel:
            arr.extend(["Flags:=", "NonModel"])

        if color is not None:
            arr.extend(["Color:=", "(%d %d %d)" % color])
        if transparency is not None:
            arr.extend(["Transparency:=", transparency])
        if material is not None:
            arr.extend(["MaterialName:=", material])
        if solve_inside is not None:
            arr.extend(["SolveInside:=", solve_inside])
        return arr

    def _selections_array(self, *names):
        return ["NAME:Selections", "Selections:=", ",".join(names)]

    def draw_box_corner(self, pos, size, **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        name = self._modeler.CreateBox(
            ["NAME:BoxParameters",
             "XPosition:=", pos[0], # '('+str(pos[0])+')'+self.unit,
             "YPosition:=", pos[1], # '('+str(pos[1])+')'+self.unit,
             "ZPosition:=", pos[2], # '('+str(pos[2])+')'+self.unit,
             "XSize:=", size[0], # '('+str(size[0])+')'+self.unit,
             "YSize:=", size[1], # '('+str(size[1])+')'+self.unit,
             "ZSize:=", size[2]], #'('+str(size[2])+')'+self.unit],
            self._attributes_array(**kwargs)
        )
        if name != kwargs['name']:
            msg = ('Warning: \'%s\' already exists, '
                   'new box named \'%s\'' % (kwargs['name'], name))
            print(msg)
        return Box(name, self, pos, size)
    
    
#    def draw_cylinder(self, pos, size, **kwargs):
#        pos = parse_entry(pos)
#        size = parse_entry(size)
#        name = self._modeler.CreateCylinder(
#            ["NAME:CylinderParameters",
#             "XCenter:=", pos[0], # '('+str(pos[0])+')'+self.unit,
#             "YCenter:=", pos[1], # '('+str(pos[1])+')'+self.unit,
#             "ZCenter:=", pos[2], # '('+str(pos[2])+')'+self.unit,
#             "Radius:=", size[0], # '('+str(size[0])+')'+self.unit,
#             "Height:=", size[1],
#             "WhichAxis:=", "Z"], # '('+str(size[1])+')'+self.unit,
#            self._attributes_array(**kwargs)
#        )
#        if name != kwargs['name']:
#            msg = ('Warning: \'%s\' already exists, '
#                   'new box named \'%s\'' % (kwargs['name'], name))
#            print(msg)
#        return Box(name, self, pos, size)

    def draw_box_center(self, pos, size, **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        corner_pos = [var(p) - var(s)/2 for p, s in zip(pos, size)]
#        print(corner_pos)
        return self.draw_box_corner(corner_pos, size, **kwargs)

    def draw_polyline(self, points, closed=True, **kwargs):
        points = parse_entry(points)
#        print(points)
        pointsStr = ["NAME:PolylinePoints"]
        indexsStr = ["NAME:PolylineSegments"]
        for ii, point in enumerate(points):
            pointsStr.append(["NAME:PLPoint", "X:=", point[0], "Y:=", point[1], "Z:=", point[2]])
            indexsStr.append(["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", ii, "NoOfPoints:=", 2])
        if closed:
            pointsStr.append(["NAME:PLPoint", "X:=", points[0][0], "Y:=", points[0][1], "Z:=", points[0][2]])
            params_closed = ["IsPolylineCovered:=", True, "IsPolylineClosed:=", True]
        else:
            indexsStr = indexsStr[:-1]
            params_closed = ["IsPolylineCovered:=", True, "IsPolylineClosed:=", False]

        name = self._modeler.CreatePolyline(
            ["NAME:PolylineParameters",
            *params_closed,
            pointsStr,
            indexsStr],
            self._attributes_array(**kwargs)
        )

        if closed:
            return Polyline(name, self, points)
        else:
            return OpenPolyline(name, self, points)

    def draw_rect_corner(self, pos, size=[0,0,0], **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        assert ('0' in size or 0 in size)
        axis = "XYZ"[size.index(0)]
        w_idx, h_idx = {'X': (1, 2),
                        'Y': (2, 0),
                        'Z': (0, 1)}[axis]
        name = self._modeler.CreateRectangle(
            ["NAME:RectangleParameters",
             "XStart:=", pos[0],
             "YStart:=", pos[1],
             "ZStart:=", pos[2],
             "Width:=", size[w_idx],
             "Height:=", size[h_idx],
             "WhichAxis:=", axis],
            self._attributes_array(**kwargs)
        )
        return Rect(name, self, pos, size)

    def draw_rect_center(self, pos, size=[0,0,0], **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        corner_pos = [var(p) - var(s)/2 for p, s in zip(pos, size)]
        return self.draw_rect_corner(corner_pos, size, **kwargs)

    def draw_cylinder(self, pos, radius, height, axis, **kwargs):
        assert axis in "XYZ"
        return self._modeler.CreateCylinder(
            ["NAME:CylinderParameters",
             "XCenter:=", pos[0],
             "YCenter:=", pos[1],
             "ZCenter:=", pos[2],
             "Radius:=", radius,
             "Height:=", height,
             "WhichAxis:=", axis,
             "NumSides:=", 0],
            self._attributes_array(**kwargs))

    def draw_disk(self, pos, radius, axis, **kwargs):
            assert axis in "XYZ"
            return self._modeler.CreateEllipse(
                ["NAME:EllipsdeParameters",
                 "XCenter:=", pos[0],
                 "YCenter:=", pos[1],
                 "ZCenter:=", pos[2],
                 "MajRadius:=", radius,
                 "Ratio:=", 1,
                 "WhichAxis:=", axis],
                self._attributes_array(**kwargs))

    def draw_cylinder_center(self, pos, radius, height, axis, **kwargs):
        assert axis in "XYZ"
        axis_idx = ["X", "Y", "Z"].index(axis)
        edge_pos = copy(pos)
        edge_pos[axis_idx] = var(pos[axis_idx]) - var(height)/2
        return self.draw_cylinder(edge_pos, radius, height, axis, **kwargs)

    def draw_wirebond(self, pos, ori, width, height='0.1mm', **kwargs): #ori should be normed
        pos, ori, width, heigth = parse_entry((pos, ori, width, height))
        xpad = pos[0]-width/2.*ori[1]
        ypad = pos[1]+width/2.*ori[0]
        xdir = ori[1]
        ydir = -ori[0]
        name = self._modeler.CreateBondwire(["NAME:BondwireParameters",
                                            "WireType:=", "Low",
                                            "WireDiameter:=", "0.02mm",
                                            "NumSides:=", 6,
                                            "XPadPos:=", xpad,
                                            "YPadPos:=", ypad,
                                            "ZPadPos:=", "0mm",
                                            "XDir:=", xdir,
                                            "YDir:=", ydir,
                                            "ZDir:=", 0,
                                            "Distance:=", width,
                                            "h1:=", height,
                                            "h2:=", "0mm",
                                            "alpha:=", "80deg",
                                            "beta:=", "80deg",
                                            "WhichAxis:=", "Z"],
                                            self._attributes_array(**kwargs))

        return name # TODO create Wirebond class
    
    def connect_faces(self, name1, name2):
        name = self._modeler.Connect(["NAME:Selections",
                               "Selections:=", ','.join([name1, name2])])
#        print(name) #None
        return name

    def delete(self, name):
        self._modeler.Delete(["NAME:Selections", "Selections:=", name])
        
    def unite(self, names, name=None, keep_originals=False):
        self._modeler.Unite(
            self._selections_array(*names),
            ["NAME:UniteParameters", "KeepOriginals:=", keep_originals]
        )
        if name is None:
            return names[0]
        else:
#            return names[0].rename(name)
            return self.rename_obj(names[0], name)

    def intersect(self, names, keep_originals=False):
        self._modeler.Intersect(
            self._selections_array(*names),
            ["NAME:IntersectParameters", "KeepOriginals:=", keep_originals]
        )
        return names[0]

    def subtract(self, blank_name, tool_names, keep_originals=False):
        selection_array= ["NAME:Selections",
                          "Blank Parts:=", blank_name,
                          "Tool Parts:=", ",".join(tool_names)]
        self._modeler.Subtract(
            selection_array,
            ["NAME:UniteParameters", "KeepOriginals:=", keep_originals]
        )
        return blank_name

    def translate(self, name, vector):
        self._modeler.Move(
            self._selections_array(name),
            ["NAME:TranslateParameters",
        		"TranslateVectorX:="	, vector[0],
        		"TranslateVectorY:="	, vector[1],
        		"TranslateVectorZ:="	, vector[2]]
        )
    
    def relative_CS(self, position, name):
        self._modeler.CreateRelativeCS(
	[
		"NAME:RelativeCSParameters",
		"Mode:="		, "Axis/Position",
		"OriginX:="		, position[0],
		"OriginY:="		, position[1],
		"OriginZ:="		, position[2],
		"XAxisXvec:="		, "1mm",
		"XAxisYvec:="		, "0mm",
		"XAxisZvec:="		, "0mm",
		"YAxisXvec:="		, "0mm",
		"YAxisYvec:="		, "1mm",
		"YAxisZvec:="		, "0mm"
	], 
	[
		"NAME:Attributes",
		"Name:="		, name
	])
        
    def set_working_CS(self, name):
        self._modeler.SetWCS(
	[
		"NAME:SetWCS Parameter",
		"Working Coordinate System:=", name,
		"RegionDepCSOk:="	, False
	])
        
    def rotation(self, name, angle):
        self._modeler.Rotate([
		"NAME:Selections",
		"Selections:="		, name,
		"NewPartsModelFlag:="	, "Model"
	], 
	[
		"NAME:RotateParameters",
		"RotateAxis:="		, "Z",
		"RotateAngle:="		, angle+"deg"
	])


    def separate_bodies(self, name):
        self._modeler.SeparateBody(["NAME:Selections",
                                		"Selections:=", name,
                                		"NewPartsModelFlag:="	, "Model"
                                	], 
                                	[
                                		"CreateGroupsForNewObjects:=", False
                                	])

#    def make_perfect_E(self, *objects):
#        print(self._boundaries.GetBoundaries())
#        name = increment_name("PerfE", self._boundaries.GetBoundaries())
#        print(name)
#        self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", objects, "InfGroundPlane:=", False])

    def assign_perfect_E(self, obj, name='PerfE'):
#        print(obj)

        if isinstance(obj, list):
            self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", obj, "InfGroundPlane:=", False])
        else:
            name = str(obj)+'_'+name
            self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", [obj], "InfGroundPlane:=", False])
            
    def assign_perfect_E_faces(self, name):
        # this is very peculiar to cavity Si chips
        faces = list(self._modeler.GetFaceIDs(name))
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

    def assign_mesh_length(self, obj, length):#, suff = '_mesh'):
        name = str(obj)
#        print(obj)
        params = ["NAME:"+name]
        params += ["RefineInside:=", False, "Enabled:=", True]
        ######## RefineInside Should be False for planar object
        params += ["Objects:=", [str(obj)]]
        params += ["RestrictElem:=", False,
			         "RestrictLength:=",  True,
			         "MaxLength:=", length]
        self._mesh.AssignLengthOp(params)
        
    def create_object_from_face(self, name):
        faces = list(self._modeler.GetFaceIDs(name))
        faces.sort()
        face = faces[0]
        self._modeler.CreateObjectFromFaces(["NAME:Selections",
                                            "Selections:="		, name,
                                            "NewPartsModelFlag:="	, "Model"],
            ["NAME:Parameters",["NAME:BodyFromFaceToParameters","FacesToDetach:="	, [int(face)]]], 
            ["CreateGroupsForNewObjects:=", False])
        return name+'_ObjectFromFace1'

    def _fillet(self, radius, vertex_index, obj):
        vertices = self._modeler.GetVertexIDsFromObject(obj)
        if isinstance(vertex_index, list):
            to_fillet = [int(vertices[v]) for v in vertex_index]
        else:
            to_fillet = [int(vertices[vertex_index])]
#        print(vertices)
#        print(radius)
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", [],
                                "Vertices:=", to_fillet,
                                "Radius:=", radius,
                                "Setback:=", "0mm"]])
            
     
    def _fillet_edges(self, radius, edge_index, obj):
        edges = self._modeler.GetEdgeIDsFromObject(obj)
        print(edges)
        if isinstance(edge_index, list):
            to_fillet = [int(edges[e]) for e in edge_index]
        else:
            to_fillet = [int(edges[edge_index])]
#        print(vertices)
#        print(radius)
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", to_fillet,
                                "Vertices:=", [],
                                "Radius:=", radius,
                                "Setback:=", "0mm"]])   
            

    def _fillets(self, radius, vertices, obj):
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", [],
                                "Vertices:=", vertices,
                                "Radius:=", radius,
                                "Setback:=", "0mm"]])

    def _sweep_along_path(self, to_sweep, path_obj):
        self.rename_obj(path_obj, str(path_obj)+'_path')
        new_name = self.rename_obj(to_sweep, path_obj)
        names = [path_obj, str(path_obj)+'_path']
        self._modeler.SweepAlongPath(self._selections_array(*names),
                                     ["NAME:PathSweepParameters",
                                		"DraftAngle:="		, "0deg",
                                		"DraftType:="		, "Round",
                                		"CheckFaceFaceIntersection:=", False,
                                		"TwistAngle:="		, "0deg"])
        return Polyline(new_name, self)
    
    
    def sweep_along_vector(self, names, vector):
        self._modeler.SweepAlongVector(self._selections_array(*names), 
                                        	["NAME:VectorSweepParameters",
                                        		"DraftAngle:="		, "0deg",
                                        		"DraftType:="		, "Round",
                                        		"CheckFaceFaceIntersection:=", False,
                                        		"SweepVectorX:="	, vector[0],
                                        		"SweepVectorY:="	, vector[1],
                                        		"SweepVectorZ:="	, vector[2]
                                        	])
                                        
                                            
    def thicken_sheet(self, sheet, thickness, bothsides=False):
        self._modeler.ThickenSheet([
                                		"NAME:Selections", "Selections:=", sheet,
                                		"NewPartsModelFlag:="	, "Model"
                                	], 
                                	[
                                		"NAME:SheetThickenParameters",
                                		"Thickness:="		, thickness,
                                		"BothSides:="		, bothsides
                                	])
                
    def mirrorZ(self, obj):
        self._modeler.Mirror([
                        		"NAME:Selections", "Selections:=", obj,
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
        return obj
    
    
    def rename_obj(self, obj, name):
        self._modeler.ChangeProperty(["NAME:AllTabs",
                                    		["NAME:Geometry3DAttributeTab",
                                    			["NAME:PropServers", str(obj)],
                                    			["NAME:ChangedProps",["NAME:Name","Value:=", str(name)]]]])
        return name


    def copy(self, obj):
        self._modeler.Copy(["NAME:Selections", "Selections:=", obj])
        new_obj = self._modeler.Paste()
        return new_obj[0]
    
    
    def duplicate_along_line(self, obj, vec):
        self._modeler.DuplicateAlongLine(["NAME:Selections","Selections:=", obj,
                                          "NewPartsModelFlag:="	, "Model"], 
                                        	["NAME:DuplicateToAlongLineParameters",
                                    		"CreateNewObjects:="	, True,
                                    		"XComponent:="		, vec[0],
                                    		"YComponent:="		, vec[1],
                                    		"ZComponent:="		, vec[2],
                                    		"NumClones:="		, "2"], 
                                        	["NAME:Options",
                                        	"DuplicateAssignments:=", False], 
                                        	["CreateGroupsForNewObjects:=", False	])

    def duplicate_along_line(self, obj, vec, n=2):
        self._modeler.DuplicateAlongLine(["NAME:Selections","Selections:=", obj,
                                          "NewPartsModelFlag:="	, "Model"], 
                                        	["NAME:DuplicateToAlongLineParameters",
                                    		"CreateNewObjects:="	, True,
                                    		"XComponent:="		, vec[0],
                                    		"YComponent:="		, vec[1],
                                    		"ZComponent:="		, vec[2],
                                    		"NumClones:="		, str(n)], 
                                        	["NAME:Options",
                                        	"DuplicateAssignments:=", False], 
                                        	["CreateGroupsForNewObjects:=", False	])    

    def _make_lumped_rlc(self, r, l, c, start, end, obj_arr, name="LumpLRC"):
        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:"+name]
        params += obj_arr
        params.append(["NAME:CurrentLine", "Start:=", start, "End:=", end])
        params += ["UseResist:=", r != 0, "Resistance:=", r,
                   "UseInduct:=", l != 0, "Inductance:=", l,
                   "UseCap:=", c != 0, "Capacitance:=", c]
        self._boundaries.AssignLumpedRLC(params)

    def _make_lumped_port(self, start, end, obj_arr, z0="50ohm", name="LumpPort"):
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
        
    def make_material(self, obj, material):
        self._modeler.ChangeProperty(["NAME:AllTabs",
                                		["NAME:Geometry3DAttributeTab",
                                			["NAME:PropServers", obj],
                                			["NAME:ChangedProps",
                                				["NAME:Material","Value:=", material]
                                			]
                                		]
                                	])
                

    def get_face_ids(self, obj):
        return self._modeler.GetFaceIDs(obj)

    def get_vertex_ids(self, obj):
        return self._modeler.GetVertexIDsFromObject(obj)

    def eval_expr(self, expr, units="mm"):
        if not isinstance(expr, str):
            return expr
        return self.parent.eval_expr(expr, units)

    def eval_var_str(self, name, unit=None):
        if not isinstance(name, VariableString):
            if unit is not None:
                return str(name)+' '+unit
            else:
                return str(name)
        return self.parent.eval_var_str(name, unit=unit)
    
    def delete_all_objects(self):
        objects = []
        for ii in range(int(self._modeler.GetNumObjects())):
            objects.append(self._modeler.GetObjectName(str(ii)))
#        print(objects)
        self._modeler.Delete(self._selections_array(*objects))
        
    def get_matched_object_name(self, name):
        return self._modeler.GetMatchedObjectName(name+'*')
    
    
class ModelEntity(str, HfssPropertyObject):
    prop_tab = "Geometry3DCmdTab"
    model_command = None
    transparency = make_float_prop("Transparent", prop_tab="Geometry3DAttributeTab", prop_server=lambda self: self)
    material = make_str_prop("Material", prop_tab="Geometry3DAttributeTab", prop_server=lambda self: self)
    coordinate_system = make_str_prop("Coordinate System")

    def __new__(self, val, *args, **kwargs):
        return str.__new__(self, val)

    def __init__(self, val, modeler):
        """
        :type val: str
        :type modeler: HfssModeler
        """
        super(ModelEntity, self).__init__()#val) #Comment out keyword to match arguments
        self.modeler = modeler
        self.prop_server = self + ":" + self.model_command + ":1"


class Box(ModelEntity):
    model_command = "CreateBox"
    position = make_float_prop("Position")
    x_size = make_float_prop("XSize")
    y_size = make_float_prop("YSize")
    z_size = make_float_prop("ZSize")
    def __init__(self, name, modeler, corner, size):
        """
        :type name: str
        :type modeler: HfssModeler
        :type corner: [(VariableString, VariableString, VariableString)]
        :param size: [(VariableString, VariableString, VariableString)]
        """
        super(Box, self).__init__(name, modeler)
        self.modeler = modeler
        self.prop_holder = modeler._modeler
        self.corner = corner
        self.size = size
        self.center = [c + s/2 for c, s in zip(corner, size)]
        faces = modeler.get_face_ids(self)
        self.z_back_face, self.z_front_face = faces[0], faces[1]
        self.y_back_face, self.y_front_face = faces[2], faces[4]
        self.x_back_face, self.x_front_face = faces[3], faces[5]

class Rect(ModelEntity):
    model_command = "CreateRectangle"
    def __init__(self, name, modeler, corner, size):
        super(Rect, self).__init__(name, modeler)
        self.corner = corner
        self.size = size
#        print(corner, size)
        self.center = [c + s/2 if s else c for c, s in zip(corner, size)]

    def make_center_line(self, axis):
        axis_idx = ["x", "y", "z"].index(axis.lower())
        start = [c for c in self.center]
        start[axis_idx] -= self.size[axis_idx]/2
        start = [self.modeler.eval_var_str(s, unit=LENGTH_UNIT) for s in start]
        end = [c for c in self.center]
        end[axis_idx] += self.size[axis_idx]/2
        end = [self.modeler.eval_var_str(s, unit=LENGTH_UNIT) for s in end]
        return start, end

    def make_rlc_boundary(self, axis, r=0, l=0, c=0, name="LumpRLC"):
        name = str(self)+'_'+name
        start, end = self.make_center_line(axis)
        self.modeler._make_lumped_rlc(r, l, c, start, end, ["Objects:=", [self]], name=name)

    def make_lumped_port(self, axis, z0="50ohm", name="LumpPort"):
        start, end = self.make_center_line(axis)
        self.modeler._make_lumped_port(start, end, ["Objects:=", [self]], z0=z0, name=name)
        
    def unite(self, list_other, name):
        union = self.modeler.unite([self] + list_other, name=name)
        return Polyline(union, self.modeler)

    def fillet(self, radius, vertex_index):
        self.modeler._fillet(radius, vertex_index, self)

    def rename(self, new_name):
        return Polyline(self.modeler.rename_obj(new_name, self), self.modeler)

class Polyline(ModelEntity): # Assume closed polyline
    model_command = "CreatePolyline"
    def __init__(self, name, modeler, points=None):
        super(Polyline, self).__init__(name, modeler)
        if points is not None:
            self.points = points
            self.n_points=len(points)
        else:
            pass
            # TODO: points = collection of points
#        axis = find_orth_axis()

# TODO: find the plane of the polyline for now, assume Z
#    def find_orth_axis():
#        X, Y, Z = (True, True, True)
#        for point in points:
#            X =
    
    def unite(self, list_other, name):
        union = self.modeler.unite(self + list_other, name=name)
        return Polyline(union, self.modeler)

    def make_center_line(self, axis): # Expects to act on a rectangle...
        #first : find center and size
        center=[0,0,0]

        for point in self.points:
            center = [center[0]+point[0]/self.n_points,
                      center[1]+point[1]/self.n_points,
                      center[2]+point[2]/self.n_points]
        size = [2*(center[0]-self.points[0][0]),
                2*(center[1]-self.points[0][1]),
                2*(center[1]-self.points[0][2])]
        axis_idx = ["x", "y", "z"].index(axis.lower())
        start = [c for c in center]
        start[axis_idx] -= size[axis_idx]/2
        start = [self.modeler.eval_var_str(s, unit=LENGTH_UNIT) for s in start]
        end = [c for c in center]
        end[axis_idx] += size[axis_idx]/2
        end = [self.modeler.eval_var_str(s, unit=LENGTH_UNIT) for s in end]
        return start, end

    def make_rlc_boundary(self, axis, r=0, l=0, c=0, name="LumpRLC"):
        name = str(self)+'_'+name
        start, end = self.make_center_line(axis)
        self.modeler._make_lumped_rlc(r, l, c, start, end, ["Objects:=", [self]], name=name)

    def fillet(self, radius, vertex_index):
        self.modeler._fillet(radius, vertex_index, self)

    def vertices(self):
        return self.modeler.get_vertex_ids(self)

    def rename(self, new_name):
        return Polyline(self.modeler.rename_obj(self, new_name), self.modeler)

class OpenPolyline(ModelEntity): # Assume closed polyline
    model_command = "CreatePolyline"
    def __init__(self, name, modeler, points=None):
        super(OpenPolyline, self).__init__(name, modeler)
        if points is not None:
            self.points = points
            self.n_points=len(points)
        else:
            pass
#        axis = find_orth_axis()

# TODO: find the plane of the polyline for now, assume Z
#    def find_orth_axis():
#        X, Y, Z = (True, True, True)
#        for point in points:
#            X =
    def vertices(self):
        return self.modeler.get_vertex_ids(self)

    def fillet(self, radius, vertex_index):
        self.modeler._fillet(radius, vertex_index, self)

    def fillets(self, radius):
        raw_list_vertices = self.modeler.get_vertex_ids(self)
        list_vertices = []
        for vertex in raw_list_vertices[1:-1]:
            list_vertices.append(int(vertex))
        if len(list_vertices)!=0:
            self.modeler._fillets(radius, list_vertices, self)
        else:
            pass

    def sweep_along_path(self, to_sweep):
        return self.modeler._sweep_along_path(to_sweep, self)

    def rename(self, new_name):
        return OpenPolyline(self.modeler.rename_obj(self, new_name), self.modeler)#, self.points)

    def copy(self, new_name):
        new_obj = OpenPolyline(self.modeler.copy(self), self.modeler)
        return new_obj.rename(new_name)


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

