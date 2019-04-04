import sys,os

class Parse:
   def __init__(self,input,debug):
      """
      Parse input topology
      """
      if not os.path.exists(input):
         raise BaseException("No input file found")
      self.options={}
      self.empty_option_dict={
                              "MolName":None,
                              "PDB":None,
                              "Ncopies":1,
                              "Placement":"random",
                              "Cutoff":0.2,
                              "Translation":[5.0,5.0,5.0],
                              "Rotation":0.5,
                              "Addon":[]
                              }
      self.empty_addon_dict={
                              "name":None,
                              "attachment":0,
                              "position":0,
                              "x0":0,
                              "k":1,
                              "N":1,
                              }
      self.input=input
      self._parse()
      self.__dict__=self.options
      
   def _parse(self):
      """
      Parse an input file
      """
      with open(self.input) as f:
         for line in f:
            if not line.lstrip().startswith("#"): # comment
               stripped_line=line.split("#")[0].strip()
               
               # Initialise an empty option dictionary with some good defaults
               if "[" in stripped_line:
                  molname=stripped_line.split()[1]
                  self.options[molname]=self.empty_option_dict.copy() # dict1=dict2 does not copy!
                  self.options[molname]["MolName"]=molname
               if ":" in stripped_line: 
                  # now process line by line
                  if "{" not in stripped_line:
                     key,value=[i.strip() for i in stripped_line.split(":")]

                     if key not in self.options[molname].keys():
                        raise BaseException("Option \"{}\" not known, please check your input file".format(key))
                     self.options[molname][key]=value                     
                  else:
                     # This is to define special lines that are given by a dictionary
                     key,value=stripped_line.split(":",1) # split on first occurence
                     if key=="Addon": # additional atoms to be added per molecule
                        addondict=self.empty_addon_dict.copy()
                        addondict_string = value.split("}",-1)[0].split("{",1)[1]
                        for pair in addondict_string.split(","):
                           addonkey,addonvalue=[i.strip() for i in pair.split(":")]
                           if addonkey not in addondict.keys():
                              raise BaseException("Option \"{}\" in Addon section of molecule {} not known, please check your input file".format(addonkey,molname))
                           addondict[addonkey]=addonvalue
                        value=addondict
                     # Since addon keyword can be used many times, this is a list
                     self.options[molname][key].append(value)                     
      self._check()
   def _check(self):
      """
      Check proper types of variables
      """
      for molname in self.options.keys():
         for key in self.options[molname].keys():
            if key in ["Ncopies"]:
               try:
                  self.options[molname][key]=int(self.options[molname][key])
               except:
                  raise BaseException("Wrong type of the variable in molecule {} section {}".format(molname,key))
            if key in ["Cutoff"]:
               try:
                  self.options[molname][key]=float(self.options[molname][key])
               except:
                  raise BaseException("Wrong type of the variable in molecule {} section {}".format(molname,key))
            if key in ["Addon"]: # test the addon part and convert variables
               for item in self.options[molname][key]: # Iterate over all attachments
                  if item is not None:
                     # attachment point
                     dtypes={"attachment":int}
                     try:
                        item["attachment"]=int(item["attachment"])
                     except:
                        raise BaseException("Wrong type of the variable in molecule {} section {}".format(molname,key))
                     # position
                     #~ try:
                        #~ print self.options[molname][key]["position"]
                        #~ self.options[molname][key]["position"]=int(self.options[molname][key]["position"])
                     #~ except:
                        #~ raise BaseException("Wrong type of the variable in molecule {} section {}".format(molname,key))
                     
   def __setitem__(self, key, item):
      self.__dict__[key] = item
   
   def __getitem__(self, key):
      return self.__dict__[key]
   
   def __repr__(self):
      return repr(self.__dict__)
   
   def __len__(self):
      return len(self.__dict__)
   
   def __delitem__(self, key):
      del self.__dict__[key]
   
   def clear(self):
      return self.__dict__.clear()
   
   def copy(self):
      return self.__dict__.copy()
   
   def has_key(self, k):
      return k in self.__dict__
   
   def update(self, *args, **kwargs):
      return self.__dict__.update(*args, **kwargs)
   
   def keys(self):
      return self.__dict__.keys()
   
   def values(self):
      return self.__dict__.values()
   
   def items(self):
      return self.__dict__.items()
   
   def pop(self, *args):
      return self.__dict__.pop(*args)
   
   def __cmp__(self, dict_):
      return self.__cmp__(self.__dict__, dict_)
   
   def __contains__(self, item):
      return item in self.__dict__
   
   def __iter__(self):
      return iter(self.__dict__)
   
   def __unicode__(self):
      return unicode(repr(self.__dict__))