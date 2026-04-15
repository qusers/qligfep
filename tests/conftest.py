"""
Pytest configuration and shared fixtures for QresFEP tests
"""

import pytest


def pytest_configure(config):
    """Register custom markers"""
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test (slow, requires Q setup)"
    )


@pytest.fixture(scope="session")
def test_data_root():
    """Session-scoped fixture for test data directory"""
    from pathlib import Path
    return Path(__file__).parent.parent / "tutorials" / "2.QresFEP_T4L"
