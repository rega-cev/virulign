#ifndef MXML_UTILS
#define MXML_UTILS

#include <vector>
#include <string>

typedef struct mxml_node_s mxml_node_t;

mxml_node_t *singleChildElement(mxml_node_t *element, 
				const std::string& tagName);

bool attributeValue(mxml_node_t *element, const std::string& attributeName,
		    std::string& result);

std::vector<mxml_node_t *> 
childElements(mxml_node_t *element, const std::string& tagName);

#endif //MXML_UTILS
