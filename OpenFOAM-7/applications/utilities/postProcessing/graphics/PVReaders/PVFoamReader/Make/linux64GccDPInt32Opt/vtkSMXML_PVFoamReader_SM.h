// Loadable modules
//
// Generated by /home/anuar/OpenFOAM/ThirdParty-7/platforms/linux64Gcc/ParaView-5.6.0/bin/vtkkwProcessXML-pv5.6
//
#ifndef __vtkSMXML_PVFoamReader_SM_h
#define __vtkSMXML_PVFoamReader_SM_h

#include <string.h>
#include <cassert>
#include <algorithm>


// From file /home/anuar/OpenFOAM/OpenFOAM-7/applications/utilities/postProcessing/graphics/PVReaders/PVFoamReader/PVFoamReader_SM.xml
static const char* const PVFoamReader_SMPVFoamReader_SMInterfaces0 =
"<ServerManagerConfiguration>\n"
"    <ProxyGroup name=\"sources\">\n"
"        <SourceProxy name=\"PVFoamReader\" class=\"vtkPVFoamReader\">\n"
"\n"
"            <!-- Send discrete time info to the animation panel -->\n"
"            <DoubleVectorProperty\n"
"                name=\"TimestepValues\"\n"
"                repeatable=\"1\"\n"
"                information_only=\"1\">\n"
"                <TimeStepsInformationHelper/>\n"
"                <Documentation>\n"
"                    Available times\n"
"                </Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <!-- File name (compulsory) -->\n"
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
"                label=\"Refresh Times\"\n"
"                animateable=\"0\">\n"
"                <Documentation>\n"
"                    Rescan for updated times and fields\n"
"                </Documentation>\n"
"            </Property>\n"
"\n"
"            <!-- Cache Mesh check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiCacheMesh\"\n"
"                command=\"SetCacheMesh\"\n"
"                label=\"Cache Mesh\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"1\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Cache the mesh in memory\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Use VTK Polyhedron check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiUseVTKPolyhedron\"\n"
"                command=\"SetUseVTKPolyhedron\"\n"
"                label=\"Use VTK Polyhedron\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Represent cells as general polyhedra instead of decomposing\n"
"                    them into simpler shapes\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Skip Zero Time check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiZeroTime\"\n"
"                command=\"SetSkipZeroTime\"\n"
"                label=\"Skip Zero Time\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Do not load the zero time directory\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Interpolate Fields check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiInterpolateVolFields\"\n"
"                command=\"SetInterpolateVolFields\"\n"
"                label=\"Interpolate Volume Fields\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"1\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Interpolate volume fields to the points\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Extrapolate Patches check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiExtrapolatePatches\"\n"
"                command=\"SetExtrapolatePatches\"\n"
"                label=\"Extrapolate Patches\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Extrapolate volume fields to non-constraint patches\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Show Patch Names check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiShowPatchNames\"\n"
"                command=\"SetShowPatchNames\"\n"
"                label=\"Show Patch Names\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Show patch names in the visualisation\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Force GUI update check box\n"
"            <IntVectorProperty\n"
"                name=\"UpdateGUI\"\n"
"                command=\"SetUpdateGUI\"\n"
"                number_of_elements=\"1\"\n"
"                is_internal=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    ???\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"                 -->\n"
"\n"
"            <!-- Include Sets check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiIncludeSets\"\n"
"                command=\"SetIncludeSets\"\n"
"                label=\"Include Sets\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <Documentation>\n"
"                    Allow selection of sets\n"
"                </Documentation>\n"
"                <BooleanDomain name=\"bool\"/>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Include Zones check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiIncludeZones\"\n"
"                command=\"SetIncludeZones\"\n"
"                label=\"Include Zones\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <Documentation>\n"
"                    Allow selection of zones\n"
"                </Documentation>\n"
"                <BooleanDomain name=\"bool\"/>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Show Groups Only check-box -->\n"
"            <IntVectorProperty\n"
"                name=\"UiShowGroupsOnly\"\n"
"                command=\"SetShowGroupsOnly\"\n"
"                label=\"Groups Only\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"0\"\n"
"                animateable=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>\n"
"                    Show only patch groups, not individual patches\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <!-- Available Parts (volume, patches, lagrangian) array -->\n"
"            <StringVectorProperty\n"
"                name=\"PartArrayStatus\"\n"
"                information_only=\"1\">\n"
"                <ArraySelectionInformationHelper attribute_name=\"Part\"/>\n"
"            </StringVectorProperty>\n"
"            <StringVectorProperty\n"
"                name=\"PartStatus\"\n"
"                label=\"Mesh Parts\"\n"
"                command=\"SetPartArrayStatus\"\n"
"                number_of_elements=\"0\"\n"
"                repeat_command=\"1\"\n"
"                number_of_elements_per_command=\"2\"\n"
"                element_types=\"2 0\"\n"
"                information_property=\"PartArrayStatus\"\n"
"                animateable=\"0\">\n"
"                <ArraySelectionDomain name=\"array_list\">\n"
"                    <RequiredProperties>\n"
"                        <Property name=\"PartArrayStatus\" function=\"ArrayList\"/>\n"
"                    </RequiredProperties>\n"
"                </ArraySelectionDomain>\n"
"                <Documentation>\n"
"                    Select the parts of the mesh to load\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <!-- Available Volume Fields array -->\n"
"            <StringVectorProperty\n"
"                name=\"VolFieldArrayStatus\"\n"
"                information_only=\"1\">\n"
"                <ArraySelectionInformationHelper attribute_name=\"VolField\"/>\n"
"            </StringVectorProperty>\n"
"            <StringVectorProperty\n"
"                name=\"VolFieldStatus\"\n"
"                label=\"Volume Fields\"\n"
"                command=\"SetVolFieldArrayStatus\"\n"
"                number_of_elements=\"0\"\n"
"                repeat_command=\"1\"\n"
"                number_of_elements_per_command=\"2\"\n"
"                element_types=\"2 0\"\n"
"                information_property=\"VolFieldArrayStatus\"\n"
"                animateable=\"0\">\n"
"                <ArraySelectionDomain name=\"array_list\">\n"
"                    <RequiredProperties>\n"
"                        <Property name=\"VolFieldArrayStatus\" function=\"ArrayList\"/>\n"
"                    </RequiredProperties>\n"
"                </ArraySelectionDomain>\n"
"                <Documentation>\n"
"                    Select volume fields to load\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <!-- Available Lagrangian fields array -->\n"
"            <StringVectorProperty\n"
"                name=\"LagrangianFieldArrayStatus\"\n"
"                information_only=\"1\">\n"
"                <ArraySelectionInformationHelper attribute_name=\"LagrangianField\"/>\n"
"            </StringVectorProperty>\n"
"            <StringVectorProperty\n"
"                name=\"LagrangianFieldStatus\"\n"
"                label=\"Lagrangian Fields\"\n"
"                command=\"SetLagrangianFieldArrayStatus\"\n"
"                number_of_elements=\"0\"\n"
"                repeat_command=\"1\"\n"
"                number_of_elements_per_command=\"2\"\n"
"                element_types=\"2 0\"\n"
"                information_property=\"LagrangianFieldArrayStatus\"\n"
"                animateable=\"0\">\n"
"                <ArraySelectionDomain name=\"array_list\">\n"
"                    <RequiredProperties>\n"
"                        <Property name=\"LagrangianFieldArrayStatus\" function=\"ArrayList\"/>\n"
"                    </RequiredProperties>\n"
"                </ArraySelectionDomain>\n"
"                <Documentation>\n"
"                    Select Lagrangain fields to load\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <!-- Available pointFields array -->\n"
"            <StringVectorProperty\n"
"                name=\"PointFieldArrayStatus\"\n"
"                information_only=\"1\">\n"
"                <ArraySelectionInformationHelper attribute_name=\"PointField\"/>\n"
"            </StringVectorProperty>\n"
"            <StringVectorProperty\n"
"                name=\"PointFieldStatus\"\n"
"                label=\"Point Fields\"\n"
"                command=\"SetPointFieldArrayStatus\"\n"
"                number_of_elements=\"0\"\n"
"                repeat_command=\"1\"\n"
"                number_of_elements_per_command=\"2\"\n"
"                element_types=\"2 0\"\n"
"                information_property=\"PointFieldArrayStatus\"\n"
"                animateable=\"0\">\n"
"                <ArraySelectionDomain name=\"array_list\">\n"
"                    <RequiredProperties>\n"
"                        <Property name=\"PointFieldArrayStatus\" function=\"ArrayList\"/>\n"
"                    </RequiredProperties>\n"
"                </ArraySelectionDomain>\n"
"                <Documentation>\n"
"                    Select point fields to load\n"
"                </Documentation>\n"
"            </StringVectorProperty>\n"
"\n"
"            <PropertyGroup label=\"Selection\">\n"
"                <Property name=\"UiIncludeSets\"/>\n"
"                <Property name=\"UiIncludeZones\"/>\n"
"                <Property name=\"UiShowGroupsOnly\"/>\n"
"                <Property name=\"PartArrayStatus\"/>\n"
"                <Property name=\"PartStatus\"/>\n"
"                <Property name=\"VolFieldArrayStatus\"/>\n"
"                <Property name=\"VolFieldStatus\"/>\n"
"                <Property name=\"LagrangianFieldArrayStatus\"/>\n"
"                <Property name=\"LagrangianFieldStatus\"/>\n"
"                <Property name=\"PointFieldArrayStatus\"/>\n"
"                <Property name=\"PointFieldStatus\"/>\n"
"            </PropertyGroup>\n"
"\n"
"            <Hints>\n"
"                <ReaderFactory\n"
"                    extensions=\"OpenFOAM\"\n"
"                    file_description=\"OpenFOAM\"/>\n"
"            </Hints>\n"
"\n"
"        </SourceProxy>\n"
"    </ProxyGroup>\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* PVFoamReader_SMPVFoamReader_SMInterfaces()
{

  const size_t len0 = strlen(PVFoamReader_SMPVFoamReader_SMInterfaces0);
  size_t len = ( 0
    + len0 );
  char* res = new char[ len + 1];
  size_t offset = 0;
  std::copy(PVFoamReader_SMPVFoamReader_SMInterfaces0, PVFoamReader_SMPVFoamReader_SMInterfaces0 + len0, res + offset); offset += len0;
  assert(offset == len);
  res[offset] = 0;
  return res;
}



#endif