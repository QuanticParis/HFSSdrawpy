# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:14:51 2019

@author: antho
"""

from KeyElement import KeyElt
from designer import Vector
from hfss import ModelEntity
import CustomElement
from hfss import parse_entry, get_active_project
from gds_modeler import GdsModeler

class PythonModeler(CustomElement.CustomElt):
    def __init__(self, interface_name,  pos=[0,0], ori=[1,0]): #"Hfss" or "Gds"
        pos1, ori1 = parse_entry((pos, ori))
        self.pos, self.ori = Vector(pos1), Vector(ori1)
        self.ports = {}
        if interface_name=="hfss":
            project = get_active_project()
            design = project.get_active_design()
            self.modeler = design.modeler
            self.modeler.set_units('mm')
            self.modeler.delete_all_objects()
            self.interface = self.modeler
            
        if interface_name=="gds":
            self.interface = GdsModeler("test")
        
    def set_units(self, units='m'):
        if (self.modeler != None):
            self.modeler.set_units('mm')
            self.interface = self.modeler

    def draw_box_corner(self, referential, pos, size, model = 'True', **kwargs):
        name = self.interface.draw_box_corner(pos, size, **kwargs)
        return ModelEntity(name,3, referential, model)
    
    def draw_box_center(self, referential, pos, size, model = 'True', **kwargs):
        name = self.interface.draw_box_center(pos, size, **kwargs)
        return (name, 3, referential, model)
  
    def draw_polyline(self, referential, pos, size, model = 'True', **kwargs):
        name = self.interface.draw_polyline(pos, size, **kwargs)
        return ModelEntity(name, 1, referential, model)
    
    def draw_rect_corner(self, referential, pos, size, model = 'True', **kwargs):
        name = self.interface.draw_rect_corner(pos, size, **kwargs)
        return ModelEntity(name, 2, referential, model)
    
    def draw_rect_center(self, referential, pos, size, model = 'True', **kwargs):
        name = self.interface.draw_rect_center(pos, size, **kwargs)
        return ModelEntity(name, 2, referential, model)
        
    def draw_cylinder(self, name, referential, pos, size, model = 'True', **kwargs):
        return ModelEntity(name, 3, referential, model, **kwargs)
    
    def draw_cylinder_center(self, name, referential, pos, size, model = 'True', **kwargs):
        return ModelEntity(name, 3, referential, pos, size, model)
    
    def draw_disk(self, name, referential, pos, size, model = 'True', **kwargs):
        return ModelEntity(name, 2, referential, model)
    
    def draw_wirebond(self, name, referential, pos, size, model = 'True', **kwargs):
        return ModelEntity(name, 2, referential, model)
    
    def connect_faces(self, name, entity1, entity2):
        assert entity1.dimension == entity2.dimension
        assert entity1.dimension == 2
        return ModelEntity(name, 3, entity1.referential, entity1.model)
        
    def delete(self, entity):
        del entity
        
    def rename_entity(self, entity, name):
        entity.rename_entity(name)
        
    def unite(self, entities, name=None, keep_originals=False):
        dim_Union = 0;
        for entity in entities:
            if entity.dimension>dim_Union:
                dim_Union = entity.dim
        if name is None:
            if dim_Union == 0:
                pass
            elif dim_Union == 1:
                union = ModelEntity(entities[0].name,1 , entities[0].referential, entities[0].model)
            elif dim_Union == 2:
                union = ModelEntity(entities[0].name,2 , entities[0].referential, entities[0].model, entities[0].boundaries)
            elif dim_Union == 3:
                union = ModelEntity(entities[0].name,3 , entities[0].referential, entities[0].model, entities[0].boundaries)
       
        else:
            if dim_Union == 0:
                pass
            elif dim_Union == 1:
                union = ModelEntity(name, 1, entities[0].referential, entities[0].model)
            elif dim_Union == 2:
                union = ModelEntity(name, 2, entities[0].referential, entities[0].model, entities[0].boundaries)
            elif dim_Union == 3:
                union = ModelEntity(name, 3, entities[0].referential, entities[0].model, entities[0].boundaries)
       
        if not(keep_originals):
            for entity in entities:
                self.delete(entity)
        return union
        
    def intersect(self, entities, keep_originals = False):
        dim_Intersection = 3;
        for entity in entities:
            if entity.dimension<dim_Intersection:
                dim_Intersection = entity.dimension
        if dim_Intersection == 0:
            pass
        elif dim_Intersection == 1:
            Intersection = ModelEntity(entities[0].name, 1, entities[0].referential, entities[0].model)
        elif dim_Intersection == 2:
            Intersection = ModelEntity(entities[0].name, 2, entities[0].referential, entities[0].model, entities[0].boundaries)
        elif dim_Intersection == 3:
            Intersection = ModelEntity(entities[0].name, 3, entities[0].referential, entities[0].model, entities[0].boundaries)
   
        #A factoriser
        if not(keep_originals):
            for entity in entities:
                self.delete(entity)
        
        return Intersection
    
    def substract(self):
        # What does this function ?
        pass

    def translate(self):
        pass
    
    def separate_bodies(self,name):
        #This looks hard
        pass
    
    def assign_perfect_E(self, entity):
        entity.boundaries.append('PerfE')
    
    def mirrorZ(self, obj):
        pass
    
    def copy(self, entity):
        if entity.dimension ==0:
            pass
        elif entity.dimension ==1:
            return ModelEntity(entity.name, entity.referential, entity.model)
        elif entity.dimension ==2:
            return ModelEntity(entity.name, entity.referential, entity.model, entity.boundaries)
        elif entity.dimension ==3:
            return ModelEntity(entity.name, entity.referential, entity.model, entity.boundaries)
        
    def duplicate_along_line(self, entity, vec):
        #Create clones
        pass
    
    def _make_lumped_rlc(self, entity):
        entity.boundaries.append('lumpedRLC')
    
    def make_material(self, entity, material):
        # Problème pour les unions de matériaux si on ajoute un attribut
        pass
    
    def delete_all_objects(self, entities):
        for entity in entities:
            self.delete(entity)
    
    
        
        
    
#    def rectangle(pos,size):
#        interface.rectangle(pos,size)
#        return ModelEntity(name)
#    
#    def chip():
#        return chip(repere,modeler)
#    

#        
#class GdsModeler():
#    def rectangle():
#        plt.plot(pos,size)
#        
class chip(PythonModeler):
    def __init__(self,repere,modeler):
            self.repere = repere
            self.modeler = modeler
#            
#chip1 = PM.chip()
#chip1.transmon(..)
#            
#            
#HfssModeler1 = HfssModeler(...)
#GdsModeler1 = GdsModeler(...)
#PM = PythonModeler(HffsModeler1)

#PM.transmon(**args)








class banque_du_labo1(KeyElt):
        
    
    def draw_JJ(self, iTrack, iGap, iTrackJ, iLength, iInduct='1nH', fillet=None):

        iTrack, iGap, iTrackJ, iLength = parse_entry((iTrack, iGap, iTrackJ, iLength))
        portOut1 = [self.coor([iLength/2,0]), self.coor_vec([1,0]), iTrack, iGap]
        portOut2 = [self.coor([-iLength/2,0]), self.coor_vec([-1,0]), iTrack, iGap]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        junction = self.connect_elt(self.name, self.name+'_1', self.name+'_2')
        pads = junction._connect_JJ(iTrackJ, iInduct=iInduct, fillet=None)
        self.trackObjects.append(pads)
        self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([iLength, iTrack+2*iGap])))


#class PythonModeler():
#    def __init__(interfaceModeler, elt_bank):
#        self.interface = interfaceModeler
#        self.elt_bank = elt_bank
#        
#PM.transmon(**args)


