#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/util/Base64.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <stack>

XERCES_CPP_NAMESPACE_USE

//Helper classes for transcoding between Unicode and local c strings
class XStr{
 public :
   XStr(const char* const toTranscode){ fUnicodeForm = XMLString::transcode(toTranscode);}
   XStr(const int number){char toTranscode[64]; sprintf(toTranscode,"%d",number);fUnicodeForm = XMLString::transcode(toTranscode);}
   XStr(const double number){char toTranscode[64]; sprintf(toTranscode,"%.16g",number);fUnicodeForm = XMLString::transcode(toTranscode);}
   
   ~XStr(){XMLString::release(&fUnicodeForm);}
   const XMLCh* unicodeForm() const{return fUnicodeForm;}
 private :
   XMLCh*   fUnicodeForm;
};
#define X(str) XStr(str).unicodeForm()

class StrX{
 public :
   StrX(const XMLCh* const toTranscode){fLocalForm = XMLString::transcode(toTranscode);}
    ~StrX(){XMLString::release(&fLocalForm);}
    const char* localForm() const{return fLocalForm;}
private :
    char*   fLocalForm;
};
inline XERCES_STD_QUALIFIER ostream& operator<<(XERCES_STD_QUALIFIER ostream& target, const StrX& toDump)
{
    target << toDump.localForm();
    return target;
}

class XMLFastNlo{
 public:
   XMLFastNlo();
   //   void PushElement(const char* const element, const char* const attname=0,const char* const attvalue=0);
   void PushElement(const char* const element, const char* const attname=0,const int attvalue=0);
   void ReplaceElement(const char* const element, const char* const attname=0,const int attvalue=0 );
   void AddValue(const char* const text);
   void AddValue(const int);
   void AddValue(const double);
   void PopElement();
   void Write(const char* const filename);
   void AddElementwithValue(const char* const element,const char* const text);
   void AddElementwithValue(const char* const element,const int number);
   void AddElementwithValue(const char* const element,const double number);
 private:
   DOMImplementation* impl;
   DOMDocument* doc;
   std::stack <DOMElement*> ElementStack;
};

XMLFastNlo::XMLFastNlo(){
   // Initialize the XML4C2 system.
   try
      {
         XMLPlatformUtils::Initialize();
      }

   catch(const XMLException& toCatch)
      {
         char *pMsg = XMLString::transcode(toCatch.getMessage());
         XERCES_STD_QUALIFIER cerr << "Error during Xerces-c Initialization.\n"
              << "  Exception message:"
              << pMsg;
         XMLString::release(&pMsg);
         exit(1);
      }
   impl =  DOMImplementationRegistry::getDOMImplementation(X("Core"));
   try
      {
         doc = impl->createDocument(
                                    0,                    // root element namespace URI.
                                    X("fastNLOtable"),         // root element name
                                    0);                   // document type object (DTD).
       }  
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
   // Initialise stack to top level
   ElementStack.push(doc->getDocumentElement());
}

// void XMLFastNlo::PushElement(const char* const element, const char* const attname,const char* const attvalue ){
//    try{
//       DOMElement *Elem = doc->createElement(X(element));
//        if(attname){
//           Elem->setAttribute(X(attname),X(attvalue));
//        }
//       (ElementStack.top())->appendChild(Elem);
//       ElementStack.push(Elem);
//    }
//    catch (const OutOfMemoryException&)
//       {
//          XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
//       }
//    catch (const XMLException& e)
//       {
//          XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
//                                    << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
//       }
// }

void XMLFastNlo::PushElement(const char* const element, const char* const attname,const int attvalue ){
   try{
      DOMElement *Elem = doc->createElement(X(element));
       if(attname){
          Elem->setAttribute(X(attname),X(attvalue));
       }
      (ElementStack.top())->appendChild(Elem);
      ElementStack.push(Elem);
   }
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
}

void XMLFastNlo::ReplaceElement(const char* const element, const char* const attname,const int attvalue ){
   try{
      DOMElement *oldElem = (DOMElement*)(ElementStack.top())->getElementsByTagName(X(element))->item(0);
      DOMElement *Elem = doc->createElement(X(element));
       if(attname){
          Elem->setAttribute(X(attname),X(attvalue));
       }
       if(oldElem){
          DOMNode *node = (ElementStack.top())->replaceChild(Elem,oldElem);
          if(node){
             node->release();
             delete node;
          }
       }else{
          (ElementStack.top())->appendChild(Elem);
       }
       ElementStack.push(Elem);
   }
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
}

void XMLFastNlo::AddValue(const char* const text){
   try{
      DOMText *DataVal = doc->createTextNode(X(text));
      (ElementStack.top())->appendChild(DataVal);
   }
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
}

void XMLFastNlo::AddValue(const int number){
   try{
      DOMText *DataVal = doc->createTextNode(X(number));
      (ElementStack.top())->appendChild(DataVal);
   }
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
}

void XMLFastNlo::AddValue(const double number){
   try{
      DOMText *DataVal = doc->createTextNode(X(number));
      (ElementStack.top())->appendChild(DataVal);
   }
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
}

void XMLFastNlo::PopElement(){
   ElementStack.pop();
}

void XMLFastNlo::Write(const char* const filename){
   try{
      DOMWriter *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
      XMLFormatTarget *myFormTarget = new LocalFileFormatTarget(X(filename));
      theSerializer->setEncoding(X("UTF-8"));
      //      theSerializer->setEncoding(X("iso-8859-1"));
      theSerializer->setFeature(X("format-pretty-print"),true);
      //   theSerializer->setNewLine(0);
      theSerializer->writeNode(myFormTarget, *doc);
      delete theSerializer;
      delete myFormTarget;
   }
   catch (const OutOfMemoryException&)
      {
         XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
      }
   catch (const XMLException& e)
      {
         XERCES_STD_QUALIFIER cerr << "An error occurred during parsing\n   Message: "
                                   << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
      }
}

void  XMLFastNlo::AddElementwithValue(const char* const element,const char* const text){
   PushElement(element);
   AddValue(text);
   PopElement();
}

void  XMLFastNlo::AddElementwithValue(const char* const element,const int number){
   PushElement(element);
   AddValue(number);
   PopElement();
}

void  XMLFastNlo::AddElementwithValue(const char* const element,const double number){
   PushElement(element);
   AddValue(number);
   PopElement();
}

