import re, struct, math

with open('/home/amrvac_bvd/amrvac-BVD/tests/mhd/icarus/icarus_2step0000.vtu', 'rb') as f:
    raw = f.read()

# Split XML header from appended binary
xml_end = raw.find(b'<AppendedData')
xml = raw[:xml_end].decode('latin-1')
# Find appended data offset marker
offset_marker = raw.find(b'_', xml_end)
appended = raw[offset_marker+1:]

# Collect all block radial ranges
blocks = []
for m in re.finditer(r'<Piece.*?NumberOfPoints="\s*(\d+)".*?NumberOfCells="\s*(\d+)"', xml):
    npts, ncells = int(m.group(1)), int(m.group(2))
    blocks.append({'npts': npts, 'ncells': ncells})

# Get first block's coordinate data (Points are first DataArray with Name="Points")
# Find the offset for Points in the first Piece
piece_starts = [m.start() for m in re.finditer(r'<Piece', xml)]
offsets = [int(m.group(1)) for m in re.finditer(r'offset="(\d+)"', xml)]

# First 3 offsets should be for Points arrays (one per block, in multi-piece)
# Actually with format="appended" each Piece has its own DataArrays with offsets
# Let me extract the first block's coordinate data
# The first Points array has the smallest offset
first_offset = min(offsets)
# Read float32 coordinates: npts * 3 floats, 4 bytes each
npts0 = blocks[0]['npts']
coord_data = struct.unpack_from(f'{npts0*3}f', appended, first_offset)

# Extract radial coords (every 3rd value starting from 0)
r_vals = sorted(set(coord_data[0::3]))
print(f'Block 1 (innermost): r=[{min(r_vals):.2f}, {max(r_vals):.2f}] R☉')
print(f'  n_cells={blocks[0]["ncells"]}, n_points={blocks[0]["npts"]}')
if len(r_vals) >= 2:
    dr_vals = [r_vals[i+1]-r_vals[i] for i in range(len(r_vals)-1)]
    print(f'  dr first: {dr_vals[0]:.3f} R☉, dr last: {dr_vals[-1]:.3f} R☉')
    if len(dr_vals) >= 2:
        q_actual = dr_vals[-1] / dr_vals[0]
        print(f'  dr_last/dr_first = {q_actual:.3f}')

# Find the outermost block (largest radial coordinate)
max_r = 0
outer_block = None
for i, off in enumerate(offsets):
    if i >= len(blocks): break
    pts = struct.unpack_from(f'{blocks[i]["npts"]*3}f', appended, off)
    rv = pts[0::3]
    rmax_block = max(rv)
    if rmax_block > max_r:
        max_r = rmax_block
        outer_block = i

if outer_block is not None:
    pts = struct.unpack_from(f'{blocks[outer_block]["npts"]*3}f', appended, offsets[outer_block])
    rv = pts[0::3]
    r_unique = sorted(set(rv))
    print(f'\nBlock {outer_block+1} (outermost): r=[{min(r_unique):.2f}, {max(r_unique):.2f}] R☉')
    if len(r_unique) >= 2:
        dr_vals = [r_unique[i+1]-r_unique[i] for i in range(len(r_unique)-1)]
        print(f'  dr first: {dr_vals[0]:.3f} R☉, dr last: {dr_vals[-1]:.3f} R☉')

print(f'\nTotal blocks: {len(blocks)}')
print(f'Radial blocks: {len(blocks)//(4*8)}  ({len(blocks)} / (4x8 theta/phi))')
