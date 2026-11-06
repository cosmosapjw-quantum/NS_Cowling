import os

def test_selected_learner_environment(run_case):
    result = run_case('environment', {})
    assert result['status'] == 'ok'
    assert result['workspace'] == os.environ.get('NS_WORKSPACE', 'solutions')
