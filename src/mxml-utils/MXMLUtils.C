#include "MXMLUtils.h"
#include "mxml/mxml.h"

#include <stdexcept>
#include <iostream>
#include <fstream>

mxml_node_t *singleChildElement(mxml_node_t *element, 
				const std::string& tagName)
{
  mxml_node_t *result = mxmlFindElement(element, element, tagName.c_str(),
					0, 0, MXML_DESCEND);

  if (result) {
    mxml_node_t *next = mxmlFindElement(result, element, tagName.c_str(),
					0, 0, MXML_NO_DESCEND);
    if (next) {
      throw std::runtime_error(std::string("Expected only one child <") 
			       + tagName
			       + "> in <" + element->value.element.name + ">");
    }
  }

  if (result && result->type != MXML_ELEMENT)
    throw std::runtime_error("Expected an XML DOM element");

  return result;
}

bool attributeValue(mxml_node_t *element, const std::string& attributeName,
		    std::string& result)
{
  const char *r = mxmlElementGetAttr(element, attributeName.c_str());

  if (r) {
    result = r;

    return true;
  } else
    return false;
}

std::vector<mxml_node_t *> 
childElements(mxml_node_t *element, const std::string& tagName)
{
  std::vector<mxml_node_t *> result;

  mxml_node_t *r = mxmlFindElement(element, element, tagName.c_str(),
				   0, 0, MXML_DESCEND);
  while (r) {
    result.push_back(r);
    r = mxmlFindElement(r, element, tagName.c_str(), 0, 0, MXML_NO_DESCEND);
  }
  
  return result;
}
