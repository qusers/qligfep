# QresFEP Test Suite

This directory contains the test suite for QresFEP.

## Structure

```
tests/
├── __init__.py                 # Package marker
├── conftest.py                 # Pytest configuration & shared fixtures
├── test_qresfep.py             # Main QresFEP test cases
└── README.md                   # This file
```

## Test Organization

Tests are organized into three tiers:

### **Tier 1: Input Validation** (`TestQresFEPValidation`)
Fast unit tests for input formats and error handling.
- Mutation format validation
- Invalid input detection
- Chain selection

### **Tier 2: File Generation & Utils** (`TestLambdaSpacing`, `TestQresFEPFileIO`)
Tests for utility functions and file operations.
- Lambda window spacing (linear, sigmoid, exponential)
- PDB file parsing
- PDB file writing

### **Tier 3: Integration Tests** (`TestQresFEPIntegration`)
Full pipeline tests using tutorial data.
- Single-topology mutation setup (A39V)
- Dual-topology mutation setup
- Multiple mutations in batch
- Different forcefield support

**Note:** Integration tests are marked with `@pytest.mark.integration` and are slower.

## Installation

### Prerequisites
```bash
pip install pytest numpy
```

### Optional
```bash
pip install pytest-cov  # For coverage reports
```

## Running Tests

### Run all tests
```bash
pytest tests/ -v
```

### Run only unit tests (fast)
```bash
pytest tests/ -v -m "not integration"
```

### Run only integration tests (slow)
```bash
pytest tests/ -v -m integration
```

### Run a specific test class
```bash
pytest tests/test_qresfep.py::TestQresFEPValidation -v
```

### Run a specific test
```bash
pytest tests/test_qresfep.py::TestLambdaSpacing::test_linear_51_windows -v
```

### Run with coverage report
```bash
pytest tests/ -v --cov=. --cov-report=html
```

## Test Data

Tests use existing tutorial data:
- **Protein PDB:** `tutorials/2.QresFEP_T4L/2LZM_prep.pdb` (prepared T4 Lysozyme)
- **Example mutations:** Defined in test cases (A39V, G25S, L153I, etc.)
- **FEP reference:** `tutorials/2.QresFEP_T4L/FEP_example/` (expected outputs)

## Next Steps

1. **Implement `Run` class fixtures:** Adapt `TestQresFEPIntegration` tests to use actual QresFEP CLI options
2. **Add PDB validation:** Extend `TestQresFEPFileIO` with full format validation
3. **Mock Q executables:** For integration tests without requiring Q installation
4. **Add parametrized tests:** Use `@pytest.mark.parametrize` for testing multiple mutations/forcefields
5. **Add coverage tracking:** Continuous integration via GitHub Actions

## Troubleshooting

### `pytest: command not found`
Install pytest: `pip install pytest`

### `Tutorial data not found`
Integration tests are skipped if tutorial files are missing. This is expected.

### `ModuleNotFoundError: No module named 'QresFEP'`
Ensure the parent directory is in `$PYTHONPATH`:
```bash
cd /workspaces/qligfep
pytest tests/ -v
```

## Contributing

When adding new features to QresFEP:
1. Add corresponding unit tests in appropriate `Test*` class
2. Mark integration tests with `@pytest.mark.integration`
3. Use fixtures from `conftest.py` for common setup
4. Keep test descriptions clear in docstrings
5. Run `pytest tests/ -v -m "not integration"` before committing
