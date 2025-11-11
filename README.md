# Smart Molecule Converter

A modern web application for converting molecular coordinate formats with professional 3D visualization capabilities.

## 🚀 Features

### Multi-Format Support
- **Z-Matrix Format** - Internal coordinates with bond lengths, angles, and dihedrals
- **Cartesian Coordinates** - Direct X, Y, Z positioning  
- **Standard XYZ Format** - Atom count with comment and coordinates
- **Mixed Format** - Combine Cartesian and Z-Matrix in the same molecule

### Smart Detection & Processing
- Automatic format detection and classification
- Cross-format referencing (Z-Matrix can reference Cartesian-defined atoms)
- Robust validation with clear error messages
- Real-time conversion with format identification

### Professional 3D Visualization
- **Advanced 3Dmol.js Integration** - High-quality molecular rendering
- **RGBA Colors** - Semi-transparent atoms with professional appearance
- **Materials Science Support** - Custom colors for Silicon, Aluminum, and Hydrogen
- **Enhanced Sizing** - Optimized atom and bond dimensions for clarity
- **Lighting Effects** - Phong shading with specular highlights
- **Atom Labels** - Element symbols and/or indices with custom positioning
- **Interactive Controls** - Smooth camera rotation, zoom, and reset functionality

### User Experience
- **Responsive Design** - Works on desktop and mobile devices
- **Real-time Feedback** - Instant conversion and visualization
- **Export Options** - Copy to clipboard and download XYZ files
- **Educational Interface** - Click-to-expand usage examples
- **Error Handling** - Clear, helpful error messages

## 🏗️ File Structure

```
/
├── index.html          # Main application interface
├── styles.css          # Styling and responsive design
├── app.js             # Core conversion logic and 3D visualization
└── README.md          # This documentation
```

## 🎯 Quick Start

1. **Open the application**
   - Open `index.html` in any modern web browser
   - No server setup required - runs entirely client-side

2. **Enter molecular coordinates**
   - Use any supported format (Z-Matrix, Cartesian, XYZ, or Mixed)
   - Click "Load Example" to see a sample formaldehyde molecule

3. **Convert and visualize**
   - Click "Convert to XYZ" to generate standard XYZ format
   - Click "View in 3D" for interactive 3D visualization
   - Toggle atom labels as needed for clarity
   - Use "Copy" or "Download" to save results

## 📖 Usage Examples

### Mixed Format (Recommended)
```
C 0.000 0.000 0.000
O 1.200 0.000 0.000
H 1 1.100 2 120.0
H 1 1.100 2 120.0 3 180.0
```

### Pure Z-Matrix
```
O
H 1 1.0100
H 1 1.0100 2 107.80
```

### Standard XYZ
```
3
Water molecule
O 0.000 0.000 0.000
H 0.757 0.000 0.000
H -0.757 0.000 0.000
```

### Cartesian Coordinates
```
O 0.000 0.000 0.000
H 0.757 0.000 0.000
H -0.757 0.000 0.000
```

## 🔬 Advanced 3D Visualization

### Professional Color Scheme
- **Hydrogen (H)**: Semi-transparent red atoms and bonds
- **Silicon (Si)**: Semi-transparent grayish gold
- **Aluminum (Al)**: Semi-transparent pink gold  
- **Other Elements**: Standard molecular colors with transparency

### Visual Effects
- **Transparency Levels**: 70-90% opacity for realistic appearance
- **Enhanced Sizing**: Hydrogen atoms 2x larger for visibility
- **Lighting**: Phong shading with specular highlights
- **Material Properties**: Realistic metallic, glass-like effects

### Customization
- **Atom Labels**: Toggle element symbols and/or indices
- **Viewer Controls**: Rotation, zoom, reset view
- **Export Options**: Copy XYZ coordinates or download files

## 🔧 Technical Details

### Format Detection Algorithm
1. Check for standard XYZ format (atom count + comment line)
2. Identify Cartesian patterns (Element + 3 numerical coordinates)
3. Detect Z-Matrix references (integer indices for bond/angle/dihedral)
4. Classify as mixed format if both patterns are present
5. Default to single format based on detected patterns

### Mixed Format Processing
1. **Parse Input** - Separate Cartesian from Z-Matrix entries
2. **Store Coordinates** - Keep all defined atomic positions
3. **Process Z-Matrix** - Convert entries using existing atoms as references
4. **Validate References** - Ensure all indices point to valid atoms
5. **Generate Output** - Combine into standard XYZ format

### Coordinate System Transformations
- **Vector Math** - Native JavaScript implementation for molecular geometry
- **Bond Lengths** - Direct distance calculations
- **Bond Angles** - Trigonometric relationships in 3D space
- **Dihedral Angles** - Complex rotations using local coordinate frames

### 3D Visualization Technology
- **3Dmol.js Library** - Professional molecular visualization
- **WebGL Rendering** - Hardware-accelerated 3D graphics
- **RGBA Color Support** - Advanced transparency and lighting effects
- **Material Science** - Custom element styling for engineering applications

## 🌐 Browser Compatibility

- **Modern Browsers** - Chrome, Firefox, Safari, Edge (latest versions)
- **JavaScript Required** - All functionality depends on JavaScript
- **WebGL Support** - Required for 3D visualization
- **No Server Required** - Runs completely in the browser

## 🧪 Key Benefits

### For Researchers
- **Flexible Input** - Use the most convenient format for each molecular fragment
- **Standard Output** - Always generate consistent XYZ format for other tools
- **Visual Verification** - See results immediately in 3D
- **Materials Science** - Support for Silicon, Aluminum, and engineering materials

### For Education
- **Progressive Learning** - Start with familiar Cartesian, add Z-Matrix complexity
- **Interactive Examples** - Click to expand detailed explanations
- **Visual Feedback** - 3D rendering helps understand molecular geometry
- **Professional Appearance** - Publication-quality molecular visualizations

### For Workflow Integration
- **Copy/Paste Ready** - Easy integration with other molecular modeling tools
- **Download Function** - Save results as standard XYZ files
- **Format Flexibility** - Handle various input formats seamlessly
- **Export Options** - Multiple ways to use converted coordinates

## 🎛️ Controls

### Main Interface
- **Load Example** - Load formaldehyde example
- **Clear** - Clear input textarea
- **Convert to XYZ** - Generate standard XYZ output
- **Center Molecule** - Option to center coordinates at origin

### 3D Viewer
- **View in 3D** - Render molecule in interactive 3D viewer
- **Reset View** - Return to default camera position
- **Show Atom Labels** - Toggle element symbols
- **Show Index Labels** - Toggle atom numbering
- **Interactive Controls** - Mouse drag to rotate, scroll to zoom

### Output
- **Copy** - Copy XYZ to clipboard
- **Download** - Save as .xyz file
- **Real-time Display** - See conversion results instantly

## 📝 Technical Notes

- **Atom Counting** - Z-Matrix references use 1-based indexing
- **Angle Units** - All angles are in degrees
- **Coordinate Precision** - Output uses 6 decimal places
- **Error Handling** - Clear messages for invalid references or formats
- **Color Format** - RGBA support with transparency levels
- **3D Performance** - Optimized for real-time molecular visualization

## 🔄 Future Enhancements

Potential improvements for future versions:
- Support for additional molecular formats (PDB, MOL, SDF)
- Batch processing of multiple molecules
- Advanced 3D visualization options (surfaces, electron density)
- Export to various molecular file formats
- Integration with molecular databases
- Animation capabilities for molecular dynamics

---

**Built for converting Z-Matrix to Cartesian coordinates and visualizing molecules in professional 3D.**