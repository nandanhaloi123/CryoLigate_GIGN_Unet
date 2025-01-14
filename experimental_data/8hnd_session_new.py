import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 17, 3, 42480])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwFOfYdVCWJhbGxTY2FsZXEDSwFHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLAUc/8AAAAAAAAH2HVQVjb2xvcnEFSwFLAH2HVQpyaWJib25UeXBlcQZLAUsAfYdVCnN0aWNrU2NhbGVxB0sBRz/wAAAAAAAAfYdVDG1tQ0lGSGVhZGVyc3EIXXEJTmFVDGFyb21hdGljTW9kZXEKSwFLAX2HVQp2ZHdEZW5zaXR5cQtLAUdAFAAAAAAAAH2HVQZoaWRkZW5xDEsBiX2HVQ1hcm9tYXRpY0NvbG9ycQ1LAU59h1UPcmliYm9uU21vb3RoaW5ncQ5LAUsAfYdVCWF1dG9jaGFpbnEPSwGIfYdVCnBkYlZlcnNpb25xEEsBSwJ9h1UIb3B0aW9uYWxxEX1xElUIb3BlbmVkQXNxE4iJSwEoVbYvbW50L2NlcGhmcy9wcm9qZWN0cy8yMDIzMTEwMTAxX0xpZ2FuZF9maXR0aW5nX3RvX0VNX21hcHMvTXlfQ29kZS9DcnlvTGlnYXRlX0dpdGh1Yi9DcnlvTGlnYXRlL0dJR04vR0lHTl9FTS9Gb3JfRWxpc2VpX05vdl8yMDI0L0NyeW9MaWdhdGVfR0lHTl9VbmV0L2V4cGVyaW1lbnRhbF9kYXRhLzhobmRfbGlnYW5kLnBkYnEUVQNQREJxFU6JdHEWfYeHc1UPbG93ZXJDYXNlQ2hhaW5zcRdLAYl9h1UJbGluZVdpZHRocRhLAUc/8AAAAAAAAH2HVQ9yZXNpZHVlTGFiZWxQb3NxGUsBSwB9h1UEbmFtZXEaSwFYDwAAADhobmRfbGlnYW5kLnBkYn2HVQ9hcm9tYXRpY0Rpc3BsYXlxG0sBiX2HVQ9yaWJib25TdGlmZm5lc3NxHEsBRz/pmZmZmZmafYdVCnBkYkhlYWRlcnNxHV1xHn1xH1gGAAAAQ1JZU1QxXXEgWEYAAABDUllTVDEgICAgMC4wMDAgICAgMC4wMDAgICAgMC4wMDAgIDkwLjAwICA5MC4wMCAgOTAuMDAgUCAxICAgICAgICAgICAxcSFhc2FVA2lkc3EiSwFLAEsAhn2HVQ5zdXJmYWNlT3BhY2l0eXEjSwFHv/AAAAAAAAB9h1UQYXJvbWF0aWNMaW5lVHlwZXEkSwFLAn2HVRRyaWJib25IaWRlc01haW5jaGFpbnElSwGIfYdVB2Rpc3BsYXlxJksBiH2HdS4='))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksBVQEgfYdVC2ZpbGxEaXNwbGF5cQNLAYl9h1UEbmFtZXEESwFYAwAAAEZZNX2HVQVjaGFpbnEFSwFYAQAAAEF9h1UOcmliYm9uRHJhd01vZGVxBksBSwJ9h1UCc3NxB0sBiYmGfYdVCG1vbGVjdWxlcQhLAUsAfYdVC3JpYmJvbkNvbG9ycQlLAU59h1UFbGFiZWxxCksBWAAAAAB9h1UKbGFiZWxDb2xvcnELSwFOfYdVCGZpbGxNb2RlcQxLAUsBfYdVBWlzSGV0cQ1LAYl9h1ULbGFiZWxPZmZzZXRxDksBTn2HVQhwb3NpdGlvbnEPXXEQTSEDSwGGcRFhVQ1yaWJib25EaXNwbGF5cRJLAYl9h1UIb3B0aW9uYWxxE31VBHNzSWRxFEsBSv////99h3Uu'))
	atomInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQdyZXNpZHVlcQJLGEsBfYdVCHZkd0NvbG9ycQNLGE59h1UEbmFtZXEESxhYAwAAAEMxOX1xBShYAwAAAEMxOF1xBksQYVgDAAAAQzEzXXEHSwFhWAMAAABDMTJdcQhLDWFYAwAAAEMxMF1xCUsAYVgDAAAAQzE3XXEKSwNhWAMAAABDMTZdcQtLD2FYAwAAAEMxNV1xDEsCYVgDAAAAQzE0XXENSw5hWAMAAABPMjRdcQ5LFmFYAwAAAE8yMl1xD0sUYVgDAAAATzIzXXEQSxVhWAMAAABPMjBdcRFLE2FYAwAAAEMwOF1xEksLYVgDAAAAQzA5XXETSwxhWAMAAABTMjFdcRRLF2FYAwAAAEMwMV1xFUsEYVgDAAAAQzAyXXEWSwVhWAMAAABDMDNdcRdLBmFYAwAAAEMwNF1xGEsHYVgDAAAAQzA1XXEZSwhhWAMAAABDMDZdcRpLCWFYAwAAAEMwN11xG0sKYVgDAAAATzExXXEcSxJhdYdVA3Zkd3EdSxiJfYdVDnN1cmZhY2VEaXNwbGF5cR5LGIl9h1UFY29sb3JxH0sYTn1xIChLAV1xIShLEksTSxRLFUsWZUsCXXEiSxdhdYdVCWlkYXRtVHlwZXEjSxiJfYdVBmFsdExvY3EkSxhVAH2HVQVsYWJlbHElSxhYAAAAAH2HVQ5zdXJmYWNlT3BhY2l0eXEmSxhHv/AAAAAAAAB9h1UHZWxlbWVudHEnSxhLBn1xKChLCF1xKShLEksTSxRLFUsWZUsQXXEqSxdhdYdVCmxhYmVsQ29sb3JxK0sYTn2HVQxzdXJmYWNlQ29sb3JxLEsYTn2HVQ9zdXJmYWNlQ2F0ZWdvcnlxLUsYWAQAAABtYWlufYdVBnJhZGl1c3EuSxhHP/4UeuAAAAB9cS8oRz/5wo9gAAAAXXEwKEsASwJLDksQZUc/9rhR4AAAAF1xMShLEksUSxVLFmVHP/dcKQAAAABdcTJLE2FHP/xR64AAAABdcTNLF2FHP/wo9cAAAABdcTQoSwNLD0sRZXWHVQpjb29yZEluZGV4cTVdcTZLAEsYhnE3YVULbGFiZWxPZmZzZXRxOEsYTn2HVRJtaW5pbXVtTGFiZWxSYWRpdXNxOUsYRwAAAAAAAAAAfYdVCGRyYXdNb2RlcTpLGEsCfYdVCG9wdGlvbmFscTt9cTwoVQxzZXJpYWxOdW1iZXJxPYiIXXE+SwFLGIZxP2GHVQdiZmFjdG9ycUCIiUsYR0BRHCkAAAAAfYeHVQlvY2N1cGFuY3lxQYiJSxhHP/AAAAAAAAB9h4d1VQdkaXNwbGF5cUJLGIh9h3Uu'))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECSxtOfYdVBWF0b21zcQNdcQQoXXEFKEsCSwdlXXEGKEsCSw5lXXEHKEsCSxRlXXEIKEsDSw9lXXEJKEsDSxBlXXEKKEsESwplXXELKEsESxBlXXEMKEsESxFlXXENKEsFSxFlXXEOKEsFSxJlXXEPKEsGSwdlXXEQKEsHSwhlXXERKEsHSwxlXXESKEsISwllXXETKEsJSwplXXEUKEsKSwtlXXEVKEsLSwxlXXEWKEsLSw9lXXEXKEsMSw1lXXEYKEsNSw5lXXEZKEsQSxNlXXEaKEsSSxNlXXEbKEsSSxVlXXEcKEsVSxllXXEdKEsWSxllXXEeKEsXSxllXXEfKEsYSxllZVUFbGFiZWxxIEsbWAAAAAB9h1UIaGFsZmJvbmRxIUsbiH2HVQZyYWRpdXNxIksbRz/JmZmgAAAAfYdVC2xhYmVsT2Zmc2V0cSNLG059h1UIZHJhd01vZGVxJEsbSwF9h1UIb3B0aW9uYWxxJX1VB2Rpc3BsYXlxJksbSwJ9h3Uu'))
	crdInfo = cPickle.loads(base64.b64decode('gAJ9cQFLAH1xAihLAF1xAyhHQE+9LxqfvndHQFEcOVgQYk5HQFIK8an7522HcQRHQFByj1wo9cNHQFIh2yLQ5WBHQFESn752yLSHcQVHQFDceuFHrhRHQFH68an7521HQFGH752yLQ6HcQZHQFFmVgQYk3VHQFJBaHKwIMVHQFGNDlYEGJOHcQdHQFAQ5WBBiTdHQFGq0OVgQYlHQFIpWBBiTdOHcQhHQFAnbItDlYFHQFFZWBBiTdNHQFH22RaHKwKHcQlHQFCCfvnbItFHQFE6LQ5WBBlHQFIKn752yLSHcQpHQFDFT987ZFpHQFFyfvnbItFHQFHfCj1wo9eHcQtHQFCcOVgQYk5HQFG9si0OVgRHQFGvO2RaHKyHcQxHQFBbdLxqfvpHQFGeyLQ5WBBHQFFrxqfvnbKHcQ1HQFAb52yLQ5ZHQFFhR64UeuFHQFGVkWhysCGHcQ5HQE91wo9cKPZHQFFr52yLQ5ZHQFGKHKwIMSeHcQ9HQE80GJN0vGpHQFEmVgQYk3VHQFHGuFHrhR+HcRBHQFAw1P3ztkZHQFHsKPXCj1xHQFFCTdLxqfyHcRFHQFDJaHKwIMVHQFInbItDlYFHQFE8rAgxJumHcRJHQFErU/fO2RdHQFII5WBBiTdHQFGvGp++dsmHcRNHQFFTMzMzMzNHQFJsi0OVgQZHQFFCTdLxqfyHcRRHQFEFcKPXCj1HQFJfKwIMSbpHQFEZ64UeuFKHcRVHQE+3Cj1wo9dHQFDrItDlYEJHQFJGhysCDEqHcRZHQFGMKPXCj1xHQFKp64UeuFJHQFEe+dsi0OWHcRdHQFFuhysCDEpHQFMjtkWhysFHQFF0m6XjU/iHcRhHQFHYIMSbpeNHQFLFsi0OVgRHQFGnvnbItDmHcRlHQFHwo9cKPXFHQFMfztkWhytHQFEnjU/fO2SHcRpHQFG8KPXCj1xHQFLtcKPXCj1HQFFcm6XjU/iHcRtlVQZhY3RpdmVxHEsAdXMu'))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (0, None, {}), 'display': (0, None, {}), 'id': (0, None, {}), 'vrmlString': [], 'name': (0, None, {})}
	colors = {u'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), u'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), u'Rf': ((0.8, 0, 0.34902), 1, u'default'), u'Ra': ((0, 0.490196, 0), 1, u'default'), u'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), u'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), u'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), u'Be': ((0.760784, 1, 0), 1, u'default'), u'Ba': ((0, 0.788235, 0), 1, u'default'), u'Bh': ((0.878431, 0, 0.219608), 1, u'default'), u'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), u'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), u'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), u'H': ((1, 1, 1), 1, u'default'), u'P': ((1, 0.501961, 0), 1, u'default'), u'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), u'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), u'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), u'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), u'Gd': ((0.270588, 1, 0.780392), 1, u'default'), u'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), u'Pr': ((0.85098, 1, 0.780392), 1, u'default'),
u'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), u'Pu': ((0, 0.419608, 1), 1, u'default'), u'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), u'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), u'Pa': ((0, 0.631373, 1), 1, u'default'), u'Pd': ((0, 0.411765, 0.521569), 1, u'default'), u'Cd': ((1, 0.85098, 0.560784), 1, u'default'), u'Po': ((0.670588, 0.360784, 0), 1, u'default'), u'Pm': ((0.639216, 1, 0.780392), 1, u'default'), u'Hs': ((0.901961, 0, 0.180392), 1, u'default'), u'Ho': ((0, 1, 0.611765), 1, u'default'), u'Hf': ((0.301961, 0.760784, 1), 1, u'default'), u'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), u'He': ((0.85098, 1, 1), 1, u'default'), u'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), u'Mg': ((0.541176, 1, 0), 1, u'default'), u'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default'), u'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), u'O': ((1, 0.0509804, 0.0509804), 1, u'default'), u'Mt': ((0.921569, 0, 0.14902), 1, u'default'), u'S': ((1, 1, 0.188235), 1, u'default'), u'W': ((0.129412, 0.580392, 0.839216), 1, u'default'),
u'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), u'Eu': ((0.380392, 1, 0.780392), 1, u'default'), u'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), u'Er': ((0, 0.901961, 0.458824), 1, u'default'), u'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), u'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), u'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), u'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), u'Nd': ((0.780392, 1, 0.780392), 1, u'default'), u'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), u'Np': ((0, 0.501961, 1), 1, u'default'), u'Fr': ((0.258824, 0, 0.4), 1, u'default'), u'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), u'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), u'B': ((1, 0.709804, 0.709804), 1, u'default'), u'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), u'Sr': ((0, 1, 0), 1, u'default'), u'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), u'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), u'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'), u'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'),
u'Sm': ((0.560784, 1, 0.780392), 1, u'default'), u'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), u'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), u'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), u'Sg': ((0.85098, 0, 0.270588), 1, u'default'), u'Se': ((1, 0.631373, 0), 1, u'default'), u'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), u'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), u'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), u'Ca': ((0.239216, 1, 0), 1, u'default'), u'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), u'Ce': ((1, 1, 0.780392), 1, u'default'), u'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), u'Lu': ((0, 0.670588, 0.141176), 1, u'default'), u'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), u'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), u'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), u'La': ((0.439216, 0.831373, 1), 1, u'default'), u'Li': ((0.8, 0.501961, 1), 1, u'default'), u'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), u'Tm': ((0, 0.831373, 0.321569), 1, u'default'),
u'Lr': ((0.780392, 0, 0.4), 1, u'default'), u'Th': ((0, 0.729412, 1), 1, u'default'), u'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), u'tan': ((0.823529, 0.705882, 0.54902), 1, u'default'), u'Te': ((0.831373, 0.478431, 0), 1, u'default'), u'Tb': ((0.188235, 1, 0.780392), 1, u'default'), u'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), u'Ta': ((0.301961, 0.65098, 1), 1, u'default'), u'Yb': ((0, 0.74902, 0.219608), 1, u'default'), u'Db': ((0.819608, 0, 0.309804), 1, u'default'), u'Dy': ((0.121569, 1, 0.780392), 1, u'default'), u'I': ((0.580392, 0, 0.580392), 1, u'default'), u'U': ((0, 0.560784, 1), 1, u'default'), u'Y': ((0.580392, 1, 1), 1, u'default'), u'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), u'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), u'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), u'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), u'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), u'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), u'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), u'Au': ((1, 0.819608, 0.137255), 1, u'default'),
u'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), u'In': ((0.65098, 0.458824, 0.45098), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'labelOffset': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 3, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (5, (u'tan', (0.823529, 0.705882, 0.54902, 1)), {(u'green', (0, 1, 0, 1)): [4], (u'S', (1, 1, 0.188235, 1)): [2], (u'O', (1, 0.0509804, 0.0509804, 1)): [1], (u'yellow', (1, 1, 0, 1)): [3]})
	viewerInfo = {'cameraAttrs': {'center': (66.853499980927, 72.116, 70.4685), 'fieldOfView': 48.655426633719, 'nearFar': (48.374633726482, 31.703372097385), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 70.4685}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 13.282354612621, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 0.66872357154013, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 4, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': None}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [
     {
      'cap_color': None,
      'class': 'Model_Capper_State',
      'display_style': None,
      'surface': ( 1, 0, ),
      'version': 1,
     },
     {
      'cap_color': None,
      'class': 'Model_Capper_State',
      'display_style': None,
      'surface': ( 2, 0, ),
      'version': 1,
     },
     {
      'cap_color': None,
      'class': 'Model_Capper_State',
      'display_style': None,
      'surface': ( 3, 0, ),
      'version': 1,
     },
     {
      'cap_color': None,
      'class': 'Model_Capper_State',
      'display_style': None,
      'surface': ( 4, 0, ),
      'version': 1,
     },
    ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': '8hnd_Ligand_Map_filtered_gaus_0.7.mrc',
       'path': '8hnd_Ligand_Map_filtered_gaus_0.7.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 1, 1, 0.7, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 13, 13, 9, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 13, 13, 9, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': False,
          'cap_faces': True,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': True,
          'dim_transparent_voxels': True,
          'flip_normals': False,
          'limit_voxel_count': True,
          'line_thickness': 1.0,
          'linear_interpolation': True,
          'maximum_intensity_projection': False,
          'mesh_lighting': True,
          'minimal_texture_memory': False,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': False,
          'smooth_lines': False,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': True,
          'subdivide_surface': False,
          'subdivision_levels': 1,
          'surface_smoothing': False,
          'two_sided_lighting': True,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': '56/=/|.Q*~O3Lm}IHA0\ton3!:=)L\tMe2',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1, 1, 0.7, 1, ),
          ( 1, 1, 0.7, 1, ),
          ( 1, 1, 0.7, 1, ),
         ],
        'solid_levels': [
          ( -7.9299638792870719e-05, 0, ),
          ( 1.3877934839077295, 0.99, ),
          ( 2.1910767555236816, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 1, 1, 0.7, 1, ),
         ],
        'surface_levels': [ 0.18940640604343492, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': False,
          'id': 2,
          'name': u'8hnd_Ligand_Map_filtered_gaus_0.7.mrc',
          'osl_identifier': u'#2',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 100.68032025477054,
            'rotation_axis': ( 0.7546777318747546, 0.04201073138588008, -0.6547492798474844, ),
            'translation': ( 32.41158291116735, 183.4136536453983, 8.076338882320591, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': '8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc',
       'path': '8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 1, 1, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 47, 47, 47, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 47, 47, 47, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': False,
          'cap_faces': True,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': True,
          'dim_transparent_voxels': True,
          'flip_normals': False,
          'limit_voxel_count': True,
          'line_thickness': 1.0,
          'linear_interpolation': True,
          'maximum_intensity_projection': False,
          'mesh_lighting': True,
          'minimal_texture_memory': False,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': False,
          'smooth_lines': False,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': True,
          'subdivide_surface': False,
          'subdivision_levels': 1,
          'surface_smoothing': False,
          'two_sided_lighting': True,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': "w16VISPI9f2Zq\tdvy8q6atVQ',]aly6/",
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 0.7, 1, 1, 1, ),
          ( 0.7, 1, 1, 1, ),
          ( 0.7, 1, 1, 1, ),
         ],
        'solid_levels': [
          ( -0.00047141605913639251, 0, ),
          ( 0.38408935537785294, 0.99, ),
          ( 2.258162260055542, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.7, 1, 1, 1, ),
         ],
        'surface_levels': [ 0.38408935537785294, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': False,
          'id': 3,
          'name': u'8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc',
          'osl_identifier': u'#3',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 100.92160279797513,
            'rotation_axis': ( 0.7521259245666586, 0.044413940893598965, -0.6575210988622714, ),
            'translation': ( 68.54349363467392, 92.95313714726707, 35.64950714227349, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': '8hnd_Ligand_map.mrc',
       'path': '8hnd_Ligand_map.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 0.7, 0.7, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 13, 13, 9, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 13, 13, 9, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': False,
          'cap_faces': True,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': True,
          'dim_transparent_voxels': True,
          'flip_normals': False,
          'limit_voxel_count': True,
          'line_thickness': 1.0,
          'linear_interpolation': True,
          'maximum_intensity_projection': False,
          'mesh_lighting': True,
          'minimal_texture_memory': False,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': False,
          'smooth_lines': False,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': True,
          'subdivide_surface': False,
          'subdivision_levels': 1,
          'surface_smoothing': False,
          'two_sided_lighting': True,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': '}t$2XdM}S\x0cP=Eo\x0cw3-+=K.^]mS5N+VF\\',
        'solid_brightness_factor': 1,
        'solid_colors': [
          ( 1.0, 1.0, 1.0, 1, ),
          ( 1.0, 1.0, 1.0, 1, ),
          ( 1.0, 1.0, 1.0, 1, ),
         ],
        'solid_levels': (
          ( 0.00011479089260102437, 0, ),
          ( 2.2155459235668182, 0.99, ),
          ( 3.433440685272217, 1, ),
         ),
        'solid_model': None,
        'surface_brightness_factor': 1,
        'surface_colors': [
          ( 0.7, 0.7, 0.7, 1, ),
         ],
        'surface_levels': [ 2.2155459235668182, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 1,
          'name': u'8hnd_Ligand_map.mrc',
          'osl_identifier': u'#1',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 100.68032025477054,
            'rotation_axis': ( 0.7546777318747546, 0.04201073138588008, -0.6547492798474844, ),
            'translation': ( 32.41158291116735, 183.4136536453983, 8.076338882320591, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': '8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox_resolved_with_L2_with_ReLU_k7_epoch253.mrc',
       'path': '8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox_resolved_with_L2_with_ReLU_k7_epoch253.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 0.7, 1, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 47, 47, 47, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 47, 47, 47, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': False,
          'cap_faces': True,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': True,
          'dim_transparent_voxels': True,
          'flip_normals': False,
          'limit_voxel_count': True,
          'line_thickness': 1.0,
          'linear_interpolation': True,
          'maximum_intensity_projection': False,
          'mesh_lighting': True,
          'minimal_texture_memory': False,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': False,
          'smooth_lines': False,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': True,
          'subdivide_surface': False,
          'subdivision_levels': 1,
          'surface_smoothing': False,
          'two_sided_lighting': True,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': 'a-wjK;.1ts)Oa+\tF-J\x0c1YW]fXUCc~5RF',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 0.7, 0.7, 1, 1, ),
          ( 0.7, 0.7, 1, 1, ),
          ( 0.7, 0.7, 1, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.023926800012588501, 0.99, ),
          ( 0.4430888891220093, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.7, 0.7, 1, 1, ),
         ],
        'surface_levels': [ 0.023926800012588501, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 4,
          'name': u'8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox_resolved_with_L2_with_ReLU_k7_epoch253.mrc',
          'osl_identifier': u'#4',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 101.17816366355711,
            'rotation_axis': ( 0.7557209960913929, 0.04338640023021746, -0.6534549688706148, ),
            'translation': ( 68.63882145086933, 92.96861961098361, 35.77300153436333, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
    ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_volume_dialog():
 volume_dialog_state = \
  {
   'adjust_camera': False,
   'auto_show_subregion': False,
   'box_padding': '0',
   'class': 'Volume_Dialog_State',
   'data_cache_size': '512',
   'focus_volume': '}t$2XdM}S\x0cP=Eo\x0cw3-+=K.^]mS5N+VF\\',
   'geometry': u'670x383+2562+266',
   'histogram_active_order': [ 2, 0, 1, ],
   'histogram_volumes': [ 'a-wjK;.1ts)Oa+\tF-J\x0c1YW]fXUCc~5RF', "w16VISPI9f2Zq\tdvy8q6atVQ',]aly6/", '}t$2XdM}S\x0cP=Eo\x0cw3-+=K.^]mS5N+VF\\', ],
   'immediate_update': True,
   'initial_colors': (
     ( 0.7, 0.7, 0.7, 1, ),
     ( 1, 1, 0.7, 1, ),
     ( 0.7, 1, 1, 1, ),
     ( 0.7, 0.7, 1, 1, ),
     ( 1, 0.7, 1, 1, ),
     ( 1, 0.7, 0.7, 1, ),
     ( 0.7, 1, 0.7, 1, ),
     ( 0.9, 0.75, 0.6, 1, ),
     ( 0.6, 0.75, 0.9, 1, ),
     ( 0.8, 0.8, 0.6, 1, ),
    ),
   'is_visible': True,
   'max_histograms': '3',
   'representation': 'surface',
   'selectable_subregions': False,
   'show_on_open': True,
   'show_plane': True,
   'shown_panels': [ 'Display style', 'Threshold and Color', ],
   'subregion_button': 'middle',
   'use_initial_colors': True,
   'version': 12,
   'voxel_limit_for_open': '256',
   'voxel_limit_for_plane': '256',
   'zone_radius': 2.0,
  }
 from VolumeViewer import session
 session.restore_volume_dialog_state(volume_dialog_state)

try:
  restore_volume_dialog()
except:
  reportRestoreError('Error restoring volume dialog')

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = []
	userXSections = []
	userResidueClasses = []
	residueData = [(1, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.35740674433659325, 0.6604015517481454, -0.6604015517481455), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.35740674433659325, 0.6604015517481454, 0.6604015517481455), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.25056280708573153, 0.25056280708573153, 0.9351131265310293), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (2560, 1342)
	xformMap = {0: (((0.75467773187475, 0.04201073138588, -0.65474927984748), 100.68032025477), (32.411582911167, 183.4136536454, 8.0763388823206), True)}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 53: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

