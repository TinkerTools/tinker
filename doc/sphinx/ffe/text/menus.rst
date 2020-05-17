Description of FFE Menus
========================

File Menu
---------

Force Field Explorer currently reads Tinker coordinate files (*.xyz, *.int), restart files (*.dyn), archive files (*.arc), induced dipole files (*u) and Protein Databank files (*.pdb). However, file writing has only been implemented for Tinker Cartesian coordinate files.

* Open 	
An attempt is made to open the selected file. Currently, Tinker coordinate files, Archive files, and Protein Databank files are recognized.

* Save As 	
Saves the active file under a different name. (This command is currently only supported for Cartesian coordinate Tinker files)

* Close
Close the active file.

* Close All	
Close all open files.

* Download from PubChem
Entering a valid molecule name will download the corresponding structure from the NIH PubChem database, convert atom types to the “Tiny” force field, and save the Tinker .xyz file to the NIH Downloads directory, then open it.

* Download from NCI
Entering a valid molecule name will download the corresponding structure from the NCI Cactus database, convert atom types to the “Tiny” force field, and save the Tinker .xyz file to the NIH Downloads directory, then open it.

* Download from PDB
Entering a valid 4 character accession will download the corresponding PDB file, save it to the current Force Field Explorer working directory, then open it.

* Load Restart Data
Tinker Restart files (*.dyn) contain atomic positions, velocities and two sets of accelerations (current and previous step), but no connectivity information. Loading a Restart file on top of an open Tinker coordinate file allows visualization of velocities and accelerations.

* Load Induced Data
Tinker Induced Dipole files (*.*u) contain induced dipole vectors. They are produced during a molecular dynamics routine using the Tinker “save-induced” keyword. Loading an Induced Dipole file on top of an open Tinker coordinate file allows visualization of induced dipoles.

* Exit
Exit shuts Force Field Explorer down, saving user preferences.

Selection Menu
--------------

The set of nodes currently selected are highlighted in the Molecular Tree. If the shift key his pressed while clicking on two nodes, all nodes in between are selected. Alternatively, by holding down the control key while clicking on a node, either in the Tree or Graphics window, the selection state of a node can be toggled. To use Graphics based picking, the Graphics Picking checkbox under the Picking menu must be checked.

* Select All
All open systems are selected.

* Restrict to Selections
Nodes that are not selected are made invisible.

* Highlight Selections
Selected nodes are highlighted in the Graphics window.

* Set Selection Color
Choose the color for selected atoms.

* Label Selected Atoms
Selected atoms are labeled by number and description.

* Label Selected Residues
Selected residues are labeled by their respective three letter identifiers.

* Set Label Font Size
The default label size is 12 point.

* Set Label Font Color
The default font color is white.

Display Menu
------------

We have intentionally maintained a simple, but elegant rendering mechanism. The Force Field Explorer data structure is self-rendering, meaning that each object “knows how to draw itself”. Rendering commands are “Top-Down”, meaning that commands filter down through the molecular tree, with each node responding and/or passing commands to its children. This paradigm positions Force Field Explorer to easily scale up to multi-level molecular engineering/nanotechnology modeling projects. What we want to avoid is a rendering/data architecture built upon “an array of atoms”. In our opinion, “an array of atoms” infrastructure can not achieve our molecular engineering design goals. While a graphics change is being rendered, one additional graphics command will be cached. Other commands are ignored while the 3D view catches up. Changes are applied to all selected nodes in the Molecular Tree View.

* Wireframe

* Tube
   Slow for large systems (more than a few thousand atoms).

* Spacefill

* Ball & Stick
Slow for large systems (more than a few thousand atoms).

* Invisible

* RMIN
Similar to Spacefill, but the sphere radii are set to the Lennard-Jones minima for the given atom.

* Fill
Fill rendered polygons.

* Points
Render only vertices.

* Lines
Render edges between vertices.
	
* Preferences – Radius
A scale factor applied to all rendered spheres and cylinders.

* Preferences – Wireframe Thickness
This specifies the thickness of wireframe representations in pixels.

* Preferences – Detail
Changes the number of vertices used in creating spheres and cylinders. More vertices (a higher number on the slider) create a smoother appearance, but slower rotations.

Color Menu
----------

* Monochrome
Black & White rendering.

* CPK
Atom based colors.

* Residue
Color based on amino/nucleic acid.

* Structure
Color based on secondary structure specified in PDB files. 

* Polymer
A different color is assigned to each polymer.

* Partial Charge
Coloring is based on partial charges that correspond to atom types specified in Tinker coordinate files.

* Vector Magnitude
Coloring is based on vector magnitude, where the vector may represent velocity, acceleration, induced dipole or force depending on user choice.

* User Color
Arbitrary colors can be assigned to selected atoms by choosing the “Apply User Color to Selections” menu item. The “User Color” can then be selected to invoke this scheme.

* Apply User Color to Selections
See Above.

* Set User Color
Specify the color to use when selecting the “Apply User Color to Selections” option.
 
Options Menu
------------

* Active System
Coordinate controls for the active system.
* Reset Rotation
Removes any rotations applied to the active system.
* Reset Translation
Removes any translation applied to the active system.
* Reset Rotation and Translation
Removes any rotation/translation applied to the active system.
* Rotate About Center of Mass
Specifies the active system should be rotated about its center of mass.
* Reset About Pick
Specifies the active system should be rotated about the node selected by graphics picking.

* Drag
Specify which system is controlled by mouse clicking/dragging within the Graphics window.
* Active System
Regardless of mouse position, the Active system is controlled.
* System Below Mouse
To move a system, it must be picked by clicking over it with the mouse.

* Reset Global Zoom
Reset the global zoom that has been applied to the entire scene.

* Reset Global Translation
Remove any global translation that has been applied to the entire scene.

* Reset Global Rotation
Remove any global rotation that has been applied to the entire scene.

* Reset Global View
Remove any global rotation and translation, then reset the global zoom.

* Toggle Hydrogens
Show/Hide hydrogen atoms.

* Full Screen
Make the graphics window fill the entire screen.
	
* Set Background Color
Change the color of the graphics window background. This is useful for switching to a white background before capturing images that might be printed.

Picking Menu
------------

Graphics picking can be used to select nodes anywhere in the structural hierarchy. By holding the control key while picking, the selection state of the picked node is toggled between selected/not selected. However, holding the shift button down while picking currently has no effect.

* Graphics Picking
Selecting this checkbox turns on Graphics Picking.
* Atom
Individual atoms will be picked.
* Bond
A Bond will be picked, with consecutive clicks cycling through all Bonds made by the picked atom.
* Angle
An Angle will be picked, with consecutive clicks cycling through all Angles made by the picked atom.
* Dihedral
A Dihedral will be picked, with consecutive clicks cycling through all Dihedrals made by the picked atom.
* Residue
* Polymer
* Molecule
* System
* Measure Distance
After selecting two atoms the distance between them will be reported in the Force Field Explorer console.
* Measure Angle
After selecting three atoms the angle formed between them will be reported in the Force Field Explorer console.
* Measure Dihedral
After selecting four atoms the dihedral formed between them will be reported in the Force Field Explorer console.

* Set Graphics Picking Color
Picked nodes will be highlighted in the color selected using this control.

Trajectory Menu
---------------

Tinker trajectories are stored in Archive files (*.arc), which are simply a concatenation of coordinate files (*.xyz, *.int).

* Oscillate
When the end of a trajectory is reached, it is played in the reverse direction instead of restarting at the beginning.

* Frame
Advance to an arbitrary frame of the active trajectory.

* Speed
Choose a frame rate for trajectory playback. If 20 frames per second is entered, for example, a delay of 50 msec between coordinate changes is enforced. Therefore, the actual frame rate will be slower than what is specified if rendering for the system is slow.

* Skip
Instead of displaying every frame during trajectory playback, only every nth frame will be rendered.

* Play
Play the active trajectory.

* Stop
      Stop the active trajectory.

* Step Forward
      Advance the active trajectory one frame.

* Step Back
      Rewind the active trajectory one frame.

* Reset
      Return the active trajectory to the first frame.

Simulation Menu
---------------

Tinker sampling and optimization methods are configured to start a server that listens for Force Field Explorer clients. Once a client – server connection is made, simulation information can be sent between Tinker and Force Field Explorer. Multiple Force Field Explorer clients can connect to the same Tinker simulation.

* Release Job
If Force Field Explorer is connected to a Tinker simulation, the connection will be closed.

* Connect to Local Job
Force Field Explorer will attempt to connect to a Tinker simulation executing on the current computer.

* Connect to Remote Job
Force Field Explorer will attempt to connect to a Tinker simulation executing on a remote computer.

* Set Remote Job Address
Specify the IP address of the remote computer where a Tinker simulation is executing. No domain name server is used, so the address should be of the form: XXX.XXX.XXX.XXX

Export Menu
-----------

Currently, Force Field Explorer can save screenshots uses JPEG and PNG formats. More formats will be supported when we target Java 1.5 (this release targets Java 1.4.2). In most cases, PNG gives better image quality when using Force Field Explorer than JPEG, probably because the default JPEG configuration we are using chooses small image size over preservation of detail.

* Capture Graphics
Capture the Graphics window to an image file.

* PNG
      Use the PNG image format.
      
* JPEG
      Use the JPEG image format.

> Window
Force Field Explorer components can be customized to some extent to suit user preferences. Many user configurable choices are saved between sessions.

* Reset Panes
This resets the slider between the Molecular Tree and Graphics Windows to the default position. The Graphics window is a heavy weight component (so hardware graphics acceleration can be used), while most Force Field Explorer components are light weight “Swing” components. Unfortunately, this effect prevents the user from dragging the slider and shrinking the Graphics window.

* Show Tree
            Show/Hide the Tree view.

* Show Tool Bar
            Show/Hide the Tool Bar.

* Show Global Axes
Show/Hide the Global Axes. The Global Axes are used both as a visual cue and to affect global rotations and translations.

* Show Console
   Show the command line Console.

* Java Metal Look & Feel
      This is the default appearance and should be consistent across platforms.
      
* Windows/Macintosh/Motif Look & Feel
This appearance should cause Force Field Explorer to more closely resemble the native platform.

Help Menu
---------

* About
This displays the Tinker logo and contact information.
