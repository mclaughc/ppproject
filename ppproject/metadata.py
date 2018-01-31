import copy

class Metadata:
  """The Metadata class represents a series of key/value properties associated
  with an object of some description. These objects could include audio chunks,
  images, or binary data."""

  def __init__(self):
    """Creates a new metadata class."""
    self.properties = {}

  def get_property(self, name):
    """Returns the property named name from within this object. If a property
    with this name does not exist, KeyError will be raised."""
    if (not name in self.properties):
      raise KeyError("Key '" + name + "' not found")
    return self.properties[name]

  def get_property_default(self, name, default):
    """Returns the property named name from within this object. If a property
    with this name does not exist, the value default will be returned."""
    if (not name in self.properties):
      return default
    return self.properties[name]

  def set_property(self, name, value):
    """Sets property name to value within this object."""
    self.properties[name] = value

  def clear_properties(self):
    """Removes all properties from this object."""
    self.properties.clear()

  def copy_properties(self, from_metadata, clear_properties=False):
      """Copies the properties from from_metadata to this metadata object.
      
      All properties are "deep copied", i.e. no references are maintained to the
      original object. If clear_properties is set, any existing properties set
      in this object will be removed."""
      if (clear_properties):
        self.clear_properties()

      for key, value in from_metadata.properties.items():
        self.properties[key] = copy.deepcopy(value)

  def remove_property(self, name):
    """Removes the property name from this object, if it exists."""
    if (not name in self.properties):
      return
    del self.properties[name]

  def format_string(self, template):
    """Returns a string with fields in the template string substituted with
    property names.

    Property names are represented by percent symbol in the template string. For
    example, template_for_%name%. %name% will be replaced with the property
    'name' from this object. If the property does not exist in this object,
    an empty string will be substituted instead. If the caller wishes to
    insert a percent symbol in the template string, use %%.    
    """
    out_str = ""
    search_property_name = ""
    in_property = False
    for char in template:
      if (in_property):
        if (char == '%'):
          if (len(search_property_name) > 0):
            prop_value = ""
            try:
              prop_value = str(self.get_property(search_property_name))
            except KeyError:
              pass
            out_str += prop_value
          search_property_name = ""
          in_property = False
        else:
          search_property_name += char
      else:
        if (char == '%'):
          in_property = True
        else:
          out_str += char

    # Handle unterminated property names
    if (in_property):
      out_str += '%'
      out_str += search_property_name

    return out_str

  def __str__(self):
    """Returns the string representation of this metadata object."""
    return self.properties.__str__()
