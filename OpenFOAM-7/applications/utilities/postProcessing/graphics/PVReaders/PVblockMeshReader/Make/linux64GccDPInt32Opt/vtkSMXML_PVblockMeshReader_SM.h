// Loadable modules
//
// Generated by /home/anuar/OpenFOAM/ThirdParty-7/platforms/linux64Gcc/ParaView-5.6.0/bin/vtkkwProcessXML-pv5.6
//
#ifndef __vtkSMXML_PVblockMeshReader_SM_h
#define __vtkSMXML_PVblockMeshReader_SM_h

#include <string.h>
#include <cassert>
#include <algorithm>


// From file /home/anuar/OpenFOAM/OpenFOAM-7/applications/utilities/postProcessing/graphics/PVReaders/PVblockMeshReader/PVblockMeshReader_SM.xml
static const char* const PVblockMeshReader_SMPVblockMeshReader_SMInterfaces0 =
"<ServerManagerConfiguration>\n"
"    <ProxyGroup name=\"sources\">\n"
"        <SourceProxy name=\"PVblockMeshReader\" class=\"vtkPVblockMeshReader\">\n"
"\n"
"            <!-- File name - compulsory -->\n"
"            <StringVectorProperty\n"
"                name=\"FileName\"\n"
"                command=\"SetFileName\"\n"
"                number_of_elements=\"1\"\n"
"                animateable=\"0\">\n"
"                <FileListDomain name=\"files\"/>\n"
"                <Documentation>\n"
"                    The file name\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <!-- Refresh button -->\n"
"            <Property\n"
"                name=\"UiRefresh\"\n"
"                command=\"SetRefresh\"\n"
"                label=\"Refresh\"\n"
"                animateable=\"0\">\n"
"                <Documentation>\n"
"                    Rescan for updated blockMeshDict.\n"
"                </Documentation>\n"
"            </Property>\n"
"\n"
"            <!-- Show Point Numbers check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiShowPointNumbers\"\n"
"                command=\"SetShowPointNumbers\"\n"
"                label=\"Show Point Numbers\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"1\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Show point numbers in the visualisation\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Available Blocks array -->\n"
"            <StringVectorProperty\n"
"                name=\"BlockArrayStatus\"\n"
"                information_only=\"1\">\n"
"                <ArraySelectionInformationHelper attribute_name=\"Block\"/>\n"
"            </StringVectorProperty>\n"
"            <StringVectorProperty\n"
"                name=\"BlockStatus\"\n"
"                label=\"Blocks\"\n"
"                command=\"SetBlockArrayStatus\"\n"
"                number_of_elements=\"0\"\n"
"                repeat_command=\"1\"\n"
"                number_of_elements_per_command=\"2\"\n"
"                element_types=\"2 0\"\n"
"                information_property=\"BlockArrayStatus\"\n"
"                animateable=\"0\">\n"
"                <ArraySelectionDomain name=\"array_list\">\n"
"                    <RequiredProperties>\n"
"                        <Property name=\"BlockArrayStatus\" function=\"ArrayList\"/>\n"
"                    </RequiredProperties>\n"
"                </ArraySelectionDomain>\n"
"                <Documentation>\n"
"                    Select the blocks to load\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <!-- Available CurvedEdges array -->\n"
"            <StringVectorProperty\n"
"                name=\"CurvedEdgesArrayStatus\"\n"
"                information_only=\"1\">\n"
"                <ArraySelectionInformationHelper attribute_name=\"CurvedEdges\"/>\n"
"            </StringVectorProperty>\n"
"            <StringVectorProperty\n"
"                name=\"CurvedEdgesStatus\"\n"
"                label=\"Curved Edges\"\n"
"                command=\"SetCurvedEdgesArrayStatus\"\n"
"                number_of_elements=\"0\"\n"
"                repeat_command=\"1\"\n"
"                number_of_elements_per_command=\"2\"\n"
"                element_types=\"2 0\"\n"
"                information_property=\"CurvedEdgesArrayStatus\"\n"
"                animateable=\"0\">\n"
"                <ArraySelectionDomain name=\"array_list\">\n"
"                    <RequiredProperties>\n"
"                        <Property name=\"CurvedEdgesArrayStatus\" function=\"ArrayList\"/>\n"
"                    </RequiredProperties>\n"
"                </ArraySelectionDomain>\n"
"                <Documentation>\n"
"                    Select the curved edges to load\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <PropertyGroup label=\"Selection\">\n"
"                <Property name=\"BlockArrayStatus\"/>\n"
"                <Property name=\"BlockStatus\"/>\n"
"                <Property name=\"CurvedEdgesArrayStatus\"/>\n"
"                <Property name=\"CurvedEdgesStatus\"/>\n"
"            </PropertyGroup>\n"
"\n"
"            <Hints>\n"
"                <ReaderFactory\n"
"                    extensions=\"blockMesh\"\n"
"                    file_description=\"OpenFOAM blockMesh\"/>\n"
"            </Hints>\n"
"\n"
"        </SourceProxy>\n"
"    </ProxyGroup>\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* PVblockMeshReader_SMPVblockMeshReader_SMInterfaces()
{

  const size_t len0 = strlen(PVblockMeshReader_SMPVblockMeshReader_SMInterfaces0);
  size_t len = ( 0
    + len0 );
  char* res = new char[ len + 1];
  size_t offset = 0;
  std::copy(PVblockMeshReader_SMPVblockMeshReader_SMInterfaces0, PVblockMeshReader_SMPVblockMeshReader_SMInterfaces0 + len0, res + offset); offset += len0;
  assert(offset == len);
  res[offset] = 0;
  return res;
}



#endif