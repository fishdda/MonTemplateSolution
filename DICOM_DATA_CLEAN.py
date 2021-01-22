import pydicom as dicom
import os
# import numpy as np
path = 'D:\\NBIA_HNSCC_DATA'
os.mkdir(os.path.join(path,'HNSCC_Oropharynx_Clean'))

for patient_id in os.listdir(os.path.join(path,'HNSCC_clean_oropharynx_PTVs')):
    os.mkdir(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id))
    os.mkdir(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'rtstruct')) # rtstruct
    os.mkdir(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'dicom'))    # ct dicom
    os.mkdir(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'dose'))     # dose
    os.mkdir(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'plan'))     # plan

    print('patient_id:{}'.format(patient_id))
    for dcm_name in os.listdir(os.path.join(path,'HNSCC_clean_oropharynx_PTVs',patient_id)):
        dcm = dicom.read_file(os.path.join(path,'HNSCC_clean_oropharynx_PTVs',patient_id,dcm_name),force=True)
        if dcm.Modality == 'CT':
            dcm.save_as(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'dicom',dcm_name))
        elif dcm.Modality == 'RTSTRUCT':
            dcm.save_as(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'rtstruct',dcm_name))
        elif dcm.Modality == 'RTDOSE':
            dcm.save_as(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'dose',dcm_name))
        elif dcm.Modality == 'RTPLAN':
            dcm.save_as(os.path.join(path,'HNSCC_Oropharynx_Clean',patient_id,'plan',dcm_name))


# def load_scan(path):
#     slices = []
#     for s in os.listdir(path):
#         dcm = dicom.read_file(path + '/' + s,force=True)
#         if dcm.Modality == 'CT':
#             slices.append(dcm)
#     slices.sort(key = lambda x: int(x.InstanceNumber))
#     try:
#         slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
#     except:
#         slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
        
#     for s in slices:
#         s.SliceThickness = slice_thickness
        
#     return slices

# def get_pixels_hu(scans):
#     image = np.stack([s.pixel_array for s in scans])
#     # Convert to int16 (from sometimes int16), 
#     # should be possible as values should always be low enough (<32k)
#     image = image.astype(np.int16)

#     # Set outside-of-scan pixels to 1
#     # The intercept is usually -1024, so air is approximately 0
#     image[image == -2000] = 0
    
#     # Convert to Hounsfield units (HU)
#     intercept = scans[0].RescaleIntercept
#     slope = scans[0].RescaleSlope
    
#     if slope != 1:
#         image = slope * image.astype(np.float64)
#         image = image.astype(np.int16)
        
#     image += np.int16(intercept)
    
#     return np.array(image, dtype=np.int16)

# output_path = 'D:\\NBIA_HNSCC_DATA\\HNSCC_Oropharynx\\'
# pt_path = 'D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_PTVs\\'

# for patient_id in os.listdir(pt_path):
#     print('patient_id:{}'.format(patient_id))
#     data_path = os.path.join(pt_path,patient_id)
#     patient = load_scan(data_path)
#     imgs = get_pixels_hu(patient)
#     os.mkdir(os.path.join(output_path,patient_id))
#     np.save(os.path.join(output_path,patient_id,"CT.npy"), imgs)

# lets test it
from dcmrtstruct2nii import dcmrtstruct2nii, list_rt_structs

strt_path = 'D:/NBIA_HNSCC_DATA/HNSCC_clean_oropharynx_PTVs/HNSCC-01-0001/HNSCC-01-0001_rtss.dcm'
ct_path   = 'D:/NBIA_HNSCC_DATA/HNSCC_clean_oropharynx_PTVs/HNSCC-01-0001/'
print(list_rt_structs(strt_path))

for patient_id in os.listdir('D:\\NBIA_HNSCC_DATA\\HNSCC_Oropharynx_Clean'):
    print('patient_id:{}'.format(patient_id))
    output_path = os.path.join('D:\\NBIA_HNSCC_DATA\\HNSCC_Oropharynx_Clean',patient_id,'mask')
    os.mkdir(output_path)

    dcmrtstruct2nii(os.path.join('D:\\NBIA_HNSCC_DATA\\HNSCC_Oropharynx_Clean',patient_id,'rtstruct',patient_id+'_rtss.dcm'),
                    os.path.join('D:\\NBIA_HNSCC_DATA\\HNSCC_Oropharynx_Clean',patient_id,'dicom'),
                    output_path)















# Copyright (C) 2016-2019 Matthew Jennings and Simon Biggs

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A DICOM RT Dose toolbox"""

from pymedphys._imports import matplotlib
from pymedphys._imports import numpy as np
from pymedphys._imports import plt, pydicom, scipy

from .constants import IMAGE_ORIENTATION_MAP
from .coords import xyz_axes_from_dataset
from .rtplan import get_surface_entry_point_with_fallback, require_gantries_be_zero
from .structure import pull_structure



def zyx_and_dose_from_dataset(dataset):
    x, y, z = xyz_axes_from_dataset(dataset)
    coords = (z, y, x)
    dose = dose_from_dataset(dataset)

    return coords, dose



def dose_from_dataset(ds, set_transfer_syntax_uid=True):
    """Extract the dose grid from a DICOM RT Dose file.
    """

    if set_transfer_syntax_uid:
        ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian

    dose = ds.pixel_array * ds.DoseGridScaling

    return dose



[docs]def dicom_dose_interpolate(interp_coords, dicom_dose_dataset):
    """Interpolates across a DICOM dose dataset.

    Parameters
    ----------
    interp_coords : tuple(z, y, x)
        A tuple of coordinates in DICOM order, z axis first, then y, then x
        where x, y, and z are DICOM axes.
    dose : pydicom.Dataset
        An RT DICOM Dose object
    """

    interp_z = np.array(interp_coords[0], copy=False)[:, None, None]
    interp_y = np.array(interp_coords[1], copy=False)[None, :, None]
    interp_x = np.array(interp_coords[2], copy=False)[None, None, :]

    coords, dicom_dose_dataset = zyx_and_dose_from_dataset(dicom_dose_dataset)
    interpolation = scipy.interpolate.RegularGridInterpolator(
        coords, dicom_dose_dataset
    )

    try:
        result = interpolation((interp_z, interp_y, interp_x))
    except ValueError:
        print(f"coords: {coords}")
        raise

    return result




[docs]def depth_dose(depths, dose_dataset, plan_dataset):
    """Interpolates dose for defined depths within a DICOM dose dataset.

    Since the DICOM dose dataset is in CT coordinates the corresponding
    DICOM plan is also required in order to calculate the conversion
    between CT coordinate space and depth.

    Currently, `depth_dose()` only supports a `dose_dataset` for which
    the patient orientation is HFS and that any beams in `plan_dataset`
    have gantry angle equal to 0 (head up). Depth is assumed to be
    purely in the y axis direction in DICOM coordinates.

    Parameters
    ----------
    depths : numpy.ndarray
        An array of depths to interpolate within the DICOM dose file. 0 is
        defined as the surface of the phantom using either the
        ``SurfaceEntryPoint`` parameter or a combination of
        ``SourceAxisDistance``, ``SourceToSurfaceDistance``, and
        ``IsocentrePosition``.
    dose_dataset : pydicom.dataset.Dataset
        The RT DICOM dose dataset to be interpolated
    plan_dataset : pydicom.dataset.Dataset
        The RT DICOM plan used to extract surface parameters and verify gantry
        angle 0 beams are used.
    """
    require_patient_orientation(dose_dataset, "HFS")
    require_gantries_be_zero(plan_dataset)
    depths = np.array(depths, copy=False)

    surface_entry_point = get_surface_entry_point_with_fallback(plan_dataset)
    depth_adjust = surface_entry_point.y

    y = depths + depth_adjust
    x, z = [surface_entry_point.x], [surface_entry_point.z]

    coords = (z, y, x)

    extracted_dose = np.squeeze(dicom_dose_interpolate(coords, dose_dataset))

    return extracted_dose




[docs]def profile(displacements, depth, direction, dose_dataset, plan_dataset):
    """Interpolates dose for cardinal angle horizontal profiles within a
    DICOM dose dataset.

    Since the DICOM dose dataset is in CT coordinates the corresponding
    DICOM plan is also required in order to calculate the conversion
    between CT coordinate space and depth and horizontal displacement.

    Currently, `profile()` only supports a `dose_dataset` for which
    the patient orientation is HFS and that any beams in `plan_dataset`
    have gantry angle equal to 0 (head up). Depth is assumed to be
    purely in the y axis direction in DICOM coordinates.

    Parameters
    ----------
    displacements : numpy.ndarray
        An array of displacements to interpolate within the DICOM dose
        file. 0 is defined in the DICOM z or x directions based either
        upon the ``SurfaceEntryPoint`` or the ``IsocenterPosition``
        depending on what is available within the DICOM plan file.
    depth : float
        The depth at which to interpolate within the DICOM dose file. 0 is
        defined as the surface of the phantom using either the
        ``SurfaceEntryPoint`` parameter or a combination of
        ``SourceAxisDistance``, ``SourceToSurfaceDistance``, and
        ``IsocentrePosition``.
    direction : str, one of ('inplane', 'inline', 'crossplane', 'crossline')
        Corresponds to the axis upon which to apply the displacements.
         - 'inplane' or 'inline' converts to DICOM z direction
         - 'crossplane' or 'crossline' converts to DICOM x direction
    dose_dataset : pydicom.dataset.Dataset
        The RT DICOM dose dataset to be interpolated
    plan_dataset : pydicom.dataset.Dataset
        The RT DICOM plan used to extract surface and isocentre
        parameters and verify gantry angle 0 beams are used.
    """

    require_patient_orientation(dose_dataset, "HFS")
    require_gantries_be_zero(plan_dataset)
    displacements = np.array(displacements, copy=False)

    surface_entry_point = get_surface_entry_point_with_fallback(plan_dataset)
    depth_adjust = surface_entry_point.y
    y = [depth + depth_adjust]

    if direction in ("inplane", "inline"):
        coords = (displacements + surface_entry_point.z, y, [surface_entry_point.x])
    elif direction in ("crossplane", "crossline"):
        coords = ([surface_entry_point.z], y, displacements + surface_entry_point.x)
    else:
        raise ValueError(
            "Expected direction to be equal to one of "
            "'inplane', 'inline', 'crossplane', or 'crossline'"
        )

    extracted_dose = np.squeeze(dicom_dose_interpolate(coords, dose_dataset))

    return extracted_dose



def get_dose_grid_structure_mask(
    structure_name: str,
    structure_dataset: "pydicom.Dataset",
    dose_dataset: "pydicom.Dataset",
):
    """Determines the 3D boolean mask defining whether or not a grid
    point is inside or outside of a defined structure.

    In its current implementation the dose grid and the planes upon
    which the structures are defined need to be aligned. This is due to
    the implementation only stepping through each structure plane and
    undergoing a 2D mask on the respective dose grid. In order to
    undergo a mask when the contours and dose grids do not align
    inter-slice contour interpolation would be required.

    For now, having two contours for the same structure name on a single
    slice is also not supported.

    Parameters
    ----------
    structure_name
        The name of the structure for which the mask is to be created
    structure_dataset : pydicom.Dataset
        An RT Structure DICOM object containing the respective
        structures.
    dose_dataset : pydicom.Dataset
        An RT Dose DICOM object from which the grid mask coordinates are
        determined.

    Raises
    ------
    ValueError
        If an unsupported contour is provided or the dose grid does not
        align with the structure planes.

    """
    x_dose, y_dose, z_dose = xyz_axes_from_dataset(dose_dataset)

    xx, yy = np.meshgrid(x_dose, y_dose)
    points = np.swapaxes(np.vstack([xx.ravel(), yy.ravel()]), 0, 1)

    x_structure, y_structure, z_structure = pull_structure(
        structure_name, structure_dataset
    )

    structure_z_values = []
    for item in z_structure:
        item = np.unique(item)
        if len(item) != 1:
            raise ValueError("Only one z value per contour supported")
        structure_z_values.append(item[0])

    structure_z_values = np.sort(structure_z_values)
    unique_structure_z_values = np.unique(structure_z_values)

    if np.any(structure_z_values != unique_structure_z_values):
        raise ValueError("Only one contour per slice is currently supported")

    sorted_dose_z = np.sort(z_dose)

    first_dose_index = np.where(sorted_dose_z == structure_z_values[0])[0][0]
    for i, z_val in enumerate(structure_z_values):
        dose_index = first_dose_index + i
        if structure_z_values[i] != sorted_dose_z[dose_index]:
            raise ValueError(
                "Only contours where both, there are no gaps in the "
                "z-axis of the contours, and the contour axis and dose "
                "axis, are aligned are supported."
            )

    mask_yxz = np.zeros((len(y_dose), len(x_dose), len(z_dose)), dtype=bool)

    for structure_index, z_val in enumerate(structure_z_values):
        dose_index = int(np.where(z_dose == z_val)[0])

        if z_structure[structure_index][0] != z_dose[dose_index]:
            raise ValueError("Structure and dose indices do not align")

        structure_polygon = matplotlib.path.Path(
            [
                (x_structure[structure_index][i], y_structure[structure_index][i])
                for i in range(len(x_structure[structure_index]))
            ]
        )

        # This logical "or" here is actually in place for the case where
        # there may be multiple contours on the one slice. That's not
        # going to be used at the moment however, as that case is not
        # yet supported in the logic above.
        mask_yxz[:, :, dose_index] = mask_yxz[:, :, dose_index] | (
            structure_polygon.contains_points(points).reshape(len(y_dose), len(x_dose))
        )

    mask_xyz = np.swapaxes(mask_yxz, 0, 1)
    mask_zyx = np.swapaxes(mask_xyz, 0, 2)

    return mask_zyx


def find_dose_within_structure(structure_name, structure_dataset, dose_dataset):
    dose = dose_from_dataset(dose_dataset)
    mask = get_dose_grid_structure_mask(structure_name, structure_dataset, dose_dataset)

    return dose[mask]


def create_dvh(structure, structure_dataset, dose_dataset):
    structure_dose_values = find_dose_within_structure(
        structure, structure_dataset, dose_dataset
    )
    hist = np.histogram(structure_dose_values, 100)
    freq = hist[0]
    bin_edge = hist[1]
    bin_mid = (bin_edge[1::] + bin_edge[:-1:]) / 2

    cumulative = np.cumsum(freq[::-1])
    cumulative = cumulative[::-1]
    bin_mid = np.append([0], bin_mid)

    cumulative = np.append(cumulative[0], cumulative)
    percent_cumulative = cumulative / cumulative[0] * 100

    plt.plot(bin_mid, percent_cumulative, label=structure)
    plt.title("DVH")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("Relative Volume (%)")


def require_patient_orientation(ds, patient_orientation):
    if not np.array_equal(
        ds.ImageOrientationPatient, np.array(IMAGE_ORIENTATION_MAP[patient_orientation])
    ):
        raise ValueError(
            "The supplied dataset has a patient "
            f"orientation other than {patient_orientation}."
        )

dose_path = 'D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_PTVs\\HNSCC-01-0001\\HNSCC-01-0001_rtdose.dcm'
ds = pydicom.read_file(dose_path,force=True)

def zyx_and_dose_from_dataset(dataset):
    x, y, z = xyz_axes_from_dataset(dataset)
    coords = (z, y, x)
    dose = dose_from_dataset(dataset)

    return coords, dose



def dose_from_dataset(ds, set_transfer_syntax_uid=True):
    """Extract the dose grid from a DICOM RT Dose file.
    """

    if set_transfer_syntax_uid:
        ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian

    dose = ds.pixel_array * ds.DoseGridScaling

    return dose

dose = dose_from_dataset(ds,set_transfer_syntax_uid=True)
# import os, sys, glob
# import numpy as np
# import dicom
# import matplotlib.pyplot as plt
# %matplotlib inline
# from skimage.draw import polygon

# def read_structure(structure):
#     contours = []
#         for i in range(len(structure.ROIContourSequence)):
#             contour = {}
#             contour['color'] = structure.ROIContourSequence[i].ROIDisplayColor
#             contour['number'] = structure.ROIContourSequence[i].RefdROINumber
#             contour['name'] = structure.StructureSetROISequence[i].ROIName
#             assert contour['number'] == structure.StructureSetROISequence[i].ROINumber
#             contour['contours'] = [s.ContourData for s in structure.ROIContourSequence[i].ContourSequence]
#             contours.append(contour)
#     return contours

# def get_mask(contours, slices):
#     z = [s.ImagePositionPatient[2] for s in slices]
#     pos_r = slices[0].ImagePositionPatient[1]
#     spacing_r = slices[0].PixelSpacing[1]
#     pos_c = slices[0].ImagePositionPatient[0]
#     spacing_c = slices[0].PixelSpacing[0]

#     label = np.zeros_like(image, dtype=np.uint8)
#     for con in contours:
#         num = int(con['number'])
#     for c in con['contours']:
#         nodes = np.array(c).reshape((-1, 3))
#         assert np.amax(np.abs(np.diff(nodes[:, 2]))) == 0
#         z_index = z.index(nodes[0, 2])
#         r = (nodes[:, 1] - pos_r) / spacing_r
#         c = (nodes[:, 0] - pos_c) / spacing_c
#         rr, cc = polygon(r, c)
#         label[rr, cc, z_index] = num

#         colors = tuple(np.array([con['color'] for con in contours]) / 255.0)
#     return label, colors

# train_data_path = "./DOI"
# train_patients = [os.path.join(train_data_path, name)
# for name in os.listdir(train_data_path) if os.path.isdir(os.path.join(train_data_path, name))]

# patient = train_patients[0] # Just get the first patient for demo
# for subdir, dirs, files in os.walk(patient):
# dcms = glob.glob(os.path.join(subdir, "*.dcm"))
# if len(dcms) == 1:
# structure = dicom.read_file(os.path.join(subdir, files[0]))
# contours = read_structure(structure)
# elif len(dcms) > 1:
# slices = [dicom.read_file(dcm) for dcm in dcms]
# slices.sort(key = lambda x: float(x.ImagePositionPatient[2]))
# image = np.stack([s.pixel_array for s in slices], axis=-1)
# label, colors = get_mask(contours, slices)

# # Plot to check slices, for example 50 to 59
# plt.figure(figsize=(15, 15))
# for i in range(9):
# plt.subplot(3, 3, i + 1)
# plt.imshow(image[..., i + 50], cmap="gray")
# plt.contour(label[..., i + 50], levels=[0.5, 1.5, 2.5, 3.5, 4.5], colors=colors)
# plt.axis('off')