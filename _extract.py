import json
with open('/home/amrvac_bvd/amrvac-BVD/base') as f:
    data = json.load(f)
msgs = data['messages']
for i, m in enumerate(msgs):
    if m['role'] == 'assistant':
        for c in m.get('content', []):
            if c.get('type') == 'text' and len(c.get('text','')) > 200:
                print(f'=== Message index {i}, text len {len(c["text"])} ===')
                print(c['text'][:500])
                print('...')
                print()
